import sys
import pyfaidx

fasta_file = sys.argv[1]
chip = sys.argv[2]
outfasta = sys.argv[3]
skipped_lines_log = sys.argv[4]
length_output_file = sys.argv[5]
flank_size = int(sys.argv[6]) if len(sys.argv) > 6 else 50

fasta = pyfaidx.Fasta(fasta_file)

with open(outfasta, "w") as outf, \
     open(skipped_lines_log, "w") as skipf, \
     open(length_output_file, "w") as length_outf, \
     open(chip, "r") as chp:
    
    length_outf.write("Marker\tLength\n")
    
    for line_num, line in enumerate(chp):
        line = line.strip()
        if not line:
            continue
            
        if line_num == 0:
            continue 
            
        elements = line.split("\t")
        if len(elements) < 3:
            marker_name = elements[0] if len(elements) > 0 else line
            skipf.write(f"{marker_name}\tn.a.\tn.a.\tNo_Position\n")
            continue
            
        name, chromosome, position = elements[0], elements[1], elements[2]
        
        try:
            position_int = int(position)
        except (ValueError, TypeError):
            skipf.write(f"{name}\tn.a.\tn.a.\tNo_Position\n")
            continue
        
            
        if chromosome not in fasta:
            skipf.write(f"{name}\tn.a.\tn.a.\tChr_Not_In_Reference\n")
            continue
            
        chrom_length = len(fasta[chromosome])
        start_1based = max(position_int - flank_size, 1)
        end_1based = min(position_int + flank_size, chrom_length)
        
        length_pos_start = position_int - start_1based
        sequence = str(fasta[chromosome][start_1based - 1 : end_1based]).upper()
        
        outf.write(f">{name}\n{sequence}\n")
        length_outf.write(f"{name}\t{length_pos_start}\n")

with open(outfasta, "r") as f:
    seq_count = sum(1 for line in f if line.startswith('>'))
print(f"Found {seq_count} marker sequences in {outfasta}")


