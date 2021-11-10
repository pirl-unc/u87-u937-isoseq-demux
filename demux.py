from os.path import join
from dataclasses import dataclass, field
import gzip
import yaml
from collections import Counter, defaultdict


def read_fastq(filename):
    desc = seq = quals = None

    with gzip.open(filename) as f:
        for line in f:
            line = line.decode("ascii")
            line = line.strip()
            if desc is None:
                assert line[0] in {"@", ">"}, "Unexpcted start of line '%s'" % line
                desc = line[1:]
            elif seq is None:
                seq = line
            elif line == "+":
                continue
            elif quals is None:
                quals = line
                yield (desc, seq, quals)
                desc = seq = quals = None
            else:
                assert False, "Unexpected line '%s'" % line

revcomp_translate_table = str.maketrans("ACTG", "TGAC")

def reverse_complement(s, table=revcomp_translate_table):
    return s.translate(table)[::-1]

assert reverse_complement("ACTG") == "CAGT"


@dataclass
class SampleInfo:
    barcode_to_sample : dict = field(default_factory=dict)
    barcode_to_source : dict = field(default_factory=dict)
    sample_to_barcodes : dict = field(default_factory=dict)
    sample_to_source : dict = field(default_factory=dict)
    unique_barcodes : set = field(default_factory=set)
    unique_five_prime_barcodes : set = field(default_factory=set)
    unique_three_prime_barcodes : set = field(default_factory=set)

    barcode_directions : dict = field(default_factory=dict)
    barcode_read_ends : dict = field(default_factory=dict)
    sources : set = field(default_factory=set)
    samples : set = field(default_factory=set)
    data_filenames : set = field(default_factory=set)


def parse_yaml_config(filename="config.yaml"):
    with open(filename) as f:
        d = yaml.load(f)
    print("=== Loaded %s ===" % filename)
    print(d)
    print("----------------------------")
    sources = set(d["sources"])
    source_to_sample_to_barcode_pair = d["barcodes"]
    info = SampleInfo()

    def add_barcode(barcode, source, sample, direction="forward", read_end=None):
        assert read_end in {"5p", "3p"}
        assert direction in {"forward", "reverse"}
        assert barcode not in info.unique_barcodes

        info.unique_barcodes.add(barcode)
        info.barcode_directions[barcode] = direction
        info.barcode_read_ends[barcode] = read_end
        info.barcode_to_source[barcode] = source
        info.barcode_to_sample[barcode] = sample
        if read_end == "5p":
            info.unique_five_prime_barcodes.add(barcode)
        else:
            info.unique_three_prime_barcodes.add(barcode)

    for source in sources:
        info.sources.add(source)
        sample_to_barcode_pair = source_to_sample_to_barcode_pair[source]
        for sample, barcode_pair in sample_to_barcode_pair.items():
            info.samples.add(sample)
            info.sample_to_source[sample] = source
            five_prime_barcode = barcode_pair["5p"]
            add_barcode(
                five_prime_barcode,
                source=source,
                sample=sample,
                direction="forward",
                read_end="5p")
            five_prime_barcode_rc = reverse_complement(five_prime_barcode)
            add_barcode(
                five_prime_barcode_rc,
                source=source,
                sample=sample,
                direction="reverse",
                read_end="5p")
            three_prime_barcode = barcode_pair["3p"]
            add_barcode(
                three_prime_barcode,
                source=source,
                sample=sample,
                direction="forward",
                read_end="3p")
            three_prime_barcode_rc = reverse_complement(three_prime_barcode)
            add_barcode(
                three_prime_barcode_rc,
                source=source,
                sample=sample,
                direction="reverse",
                read_end="3p")
            info.sample_to_barcodes[sample] = (
                five_prime_barcode,
                three_prime_barcode,
                five_prime_barcode_rc,
                three_prime_barcode_rc
            )
    for filename in d["data"]:
        info.data_filenames.add(filename)

    return info

def load_fastq_and_demux(
        info,
        output_dir="demux-results",
        max_trimmed_nucleotides_per_barcode_end=6):
    five_prime_barcode_lengths = {len(seq) for seq in info.unique_five_prime_barcodes}
    three_prime_barcode_lengths = {len(seq) for seq in info.unique_three_prime_barcodes}

    barcode_counts = Counter()
    sample_counts = Counter()
    if len(five_prime_barcode_lengths) != 1:
        raise ValueError("More than one barcode length")
    five_prime_length = list(five_prime_barcode_lengths)[0]
    three_prime_length = list(three_prime_barcode_lengths)[0]

    # create index for barcodes with two fewer nucleotides
    partial_barcode_to_sample = {}
    partial_barcode_to_direction = {}
    partial_barcode_blacklist = set()

    output_buffers = defaultdict(list)
    for barcode, sample in info.barcode_to_sample.items():
        direction = info.barcode_directions[barcode]
        partial_barcodes = set()
        for n in range(1, max_trimmed_nucleotides_per_barcode_end + 1):
            partial_barcodes.add(barcode[n:])
            partial_barcodes.add(barcode[:-n])

        for partial_barcode in partial_barcodes:
            if partial_barcode in partial_barcode_blacklist:
                continue
            if partial_barcode in partial_barcode_to_sample:
                del partial_barcode_to_sample[partial_barcode]
                del partial_barcode_to_direction[partial_barcode]
                partial_barcode_blacklist.add(partial_barcode)
                print("Ambiguous partial barcode: %s" % partial_barcode)
                continue
            partial_barcode_to_sample[partial_barcode] = sample
            partial_barcode_to_direction[partial_barcode] = direction
    print("Generated %d unambiguous partial barcodes, excluding %d ambiguous ones" % (
        len(partial_barcode_to_sample), len(partial_barcode_blacklist)))

    buffers = defaultdict(list)
    output_files = {}
    for sample in info.samples:
        output_files[sample] = gzip.open("%s.fa.gz" % sample, 'at', 9)

    for input_filename in info.data_filenames:
        print("Loading %s..." % input_filename)
        fastq_gen = read_fastq(input_filename)
        for (desc, seq, quals) in fastq_gen:
            sample = None
            direction = None
            barcode_5p = seq[:five_prime_length]

            sample = info.barcode_to_sample.get(barcode_5p)
            direction = info.barcode_directions.get(barcode_5p)
            if not sample:
                barcode_3p = seq[-three_prime_length:]
                sample = info.barcode_to_sample.get(barcode_3p)
                direction = info.barcode_directions.get(barcode_3p)

            if not sample:
                barcode_3p_rc = seq[:three_prime_length]

                sample = info.barcode_to_sample.get(barcode_3p_rc)
                direction = info.barcode_directions.get(barcode_3p_rc)

            if not sample:
                barcode_5p_rc = seq[-five_prime_length:]
                sample = info.barcode_to_sample.get(barcode_5p_rc)
                direction = info.barcode_directions.get(barcode_5p_rc)

            if not sample:
                for n in range(max_trimmed_nucleotides_per_barcode_end + 1):
                    for barcode in [barcode_5p, barcode_5p_rc, barcode_3p, barcode_3p_rc]:
                        barcode_trim_left = barcode[n:]
                        if barcode_trim_left in partial_barcode_to_sample:
                            sample = partial_barcode_to_sample[barcode_trim_left]
                            direction = partial_barcode_to_direction[barcode_trim_left]
                            break

                        barcode_trim_right = barcode[:-n]
                        if barcode_trim_right in partial_barcode_to_sample:
                            sample = partial_barcode_to_sample[barcode_trim_right]
                            direction = partial_barcode_to_direction[barcode_trim_right]
                            break
                    if sample:
                        break

            if sample and direction == "reverse":
                seq = reverse_complement(seq)
                quals = quals[::-1]
                barcode_5p = seq[:five_prime_length]
                barcode_3p = seq[-three_prime_length:]

            if sample:
                buffers[sample].append(">%s sample=%s direction=%s\n%s\n%s\n" % (
                    desc,
                    sample,
                    direction,
                    seq,
                    quals
                ))
                if len(buffers[sample]) > 10000:
                    output_files[sample].write("".join(buffers[sample]))
                    buffers[sample] = []
            else:
                sample = "unknown"
                direction = None
                barcode_counts[barcode_5p] += 1

            sample_counts[sample] += 1
            """
            if barcode_5p not in info.unique_five_prime_barcodes:
                print("Unknown barcode for %s" % desc)
            else:
                print("Barcode found for %s: %s" % (desc, info.barcode_to_sample[barcode_5p]))
            """
    for sample, count in sample_counts.most_common(30):
        print(sample, count)

    print("Top 5p barcodes for reads which could not be demuxed:")
    for barcode, count in barcode_counts.most_common(10):
        print(barcode, count)

    # output_path = join(output_dir, sample)

if __name__ == "__main__":
    info = parse_yaml_config()
    load_fastq_and_demux(info)
