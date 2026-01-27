import subprocess
import sys
import os

def get_chromosome_lengths(fasta_path):
    fai_path = fasta_path + ".fai" if not fasta_path.endswith(".fai") else fasta_path
    lengths = {}
    with open(fai_path) as fai:
        for line in fai:
            fields = line.strip().split('\t')
            chrom, length = fields[0], int(fields[1])
            lengths[chrom] = length
    return lengths

def fetch_flanking_sequence(fasta_path, chrom, pos, chrom_length_dict, flank):
    chrom_length = chrom_length_dict.get(chrom)
    if not chrom_length:
        raise ValueError(f"Chromosome {chrom} not found in fasta index (.fai)")
    start = max(pos - flank, 1)
    end = min(pos + flank, chrom_length)
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", fasta_path, region],
        capture_output=True,
        text=True,
        check=True
    )
    lines = result.stdout.strip().split('\n')
    seq = ''.join(lines[1:]).upper()
    return seq

# Get flank size as a command-line argument (use default 50)
try:
    fasta_file = sys.argv[1]
    chip = sys.argv[2]
    outfasta = sys.argv[3]
    skipped_lines_log = sys.argv[4]
    length_output_file = sys.argv[5]
    flank_size = int(sys.argv[6]) if len(sys.argv) > 6 else 50
except (IndexError, ValueError) as e:
    #print("Usage: python flank_seq.py <fasta> <chip> <outfasta> <outna> [flank_size]")
    sys.exit(1)

chrom_lengths = get_chromosome_lengths(fasta_file)


# Open files
with open(outfasta, "w") as outf, open(skipped_lines_log, "w") as skipf, open(length_output_file, "w") as length_outf:
    with open(chip) as chp:
        rows = [line.strip("\n") for line in chp.readlines()]
        for el in range(0, len(rows)):
            if el == 0:
                length_outf.write("Marker\tLength\n")  # header for length output
                continue
            elements = rows[el].split("\t")
            name, chromosome, position = elements[0], elements[1], elements[2]
            try:
                position_int = int(position)
            except (ValueError, TypeError):
                skipf.write(f"{rows[el]}\n")
                continue
            # Determine flank dynamically (if supporting custom flanks)
            start = max(position_int - flank_size, 1)
            length_pos_start = position_int - start
            sequence = fetch_flanking_sequence(fasta_file, chromosome, position_int, chrom_lengths, flank=flank_size)
            outf.write(f">{name}\n{sequence}\n")
            # Write pos-start with marker name to new output file
            length_outf.write(f"{name}\t{length_pos_start}\n")


print(outfasta)
# Count markers with NA and markers with derived sequences
output_lines = sum(1 for _ in open(skipped_lines_log))

with open(outfasta, "r") as f:
    seq_count = sum(1 for line in f if line.startswith('>'))
    print(f"Found {seq_count} marker sequences in {outfasta}")
    print(f"Output file {skipped_lines_log} contains {output_lines} markers with unknown position (na)")



