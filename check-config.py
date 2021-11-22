from dataclasses import dataclass, field

import yaml

@dataclass
class SampleInfo:
    barcode_to_sample  : dict = field(default_factory=dict)
    barcode_to_source : dict = field(default_factory=dict)
    sample_to_barcode_pair : dict = field(default_factory=dict)
    sample_to_source : dict = field(default_factory=dict)
    unique_barcodes : set = field(default_factory=set)
    sources : set = field(default_factory=set)
    samples : set = field(default_factory=set)

def parse_yaml(filename="config.yaml"):
    with open(filename) as f:
        d = yaml.load(f)
    print(d)
    sources = set(d["sources"])
    source_to_sample_to_barcode_pair = d["barcodes"]
    info = SampleInfo()
    for source in sources:
        info.sources.add(source)
        sample_to_barcode_pair = source_to_sample_to_barcode_pair[source]
        for sample, barcode_pair in sample_to_barcode_pair.items():
            info.samples.add(sample)
            info.sample_to_source[sample] = source
            five_prime_barcode = barcode_pair["5p"]
            assert five_prime_barcode not in info.unique_barcodes
            info.unique_barcodes.add(five_prime_barcode)

            three_prime_barcode = barcode_pair["3p"]
            assert three_prime_barcode not in info.unique_barcodes
            info.unique_barcodes.add(three_prime_barcode)

            info.barcode_to_sample[five_prime_barcode] = sample
            info.barcode_to_sample[three_prime_barcode] = sample

            info.barcode_to_source[five_prime_barcode] = source
            info.barcode_to_source[three_prime_barcode] = source

            info.sample_to_barcode_pair[sample] = (five_prime_barcode, three_prime_barcode)
    return info

if __name__ == "__main__":
    info = parse_yaml()
