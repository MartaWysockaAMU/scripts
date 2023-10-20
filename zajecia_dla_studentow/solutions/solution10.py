from Bio import SeqIO

p1, p2 = {}, {}

fastq_parser1 = SeqIO.parse("fastx_out2pe1.fastq", "fastq") 
for fastq_rec in fastq_parser1:
    p1[fastq_rec.name[:-2]] = (fastq_rec, fastq_rec.seq, fastq_rec.letter_annotations)
    
fastq_parser2 = SeqIO.parse("fastx_out2pe2.fastq", "fastq") 
for fastq_rec in fastq_parser2:
    p2[fastq_rec.name[:-2]] = (fastq_rec, fastq_rec.seq, fastq_rec.letter_annotations)

handle_paired1 = open("paired_p1.fastq", "w")
handle_paired2 = open("paired_p2.fastq", "w")
handle_single1 = open("single_p1.fastq", "w")
handle_single2 = open("single_p2.fastq", "w")
for key, value in p1.items():
    if key in p2.keys():
        SeqIO.write( (p1[key])[0], handle_paired1, "fastq")
        SeqIO.write( (p2[key])[0], handle_paired2, "fastq")
        # Ma pare. Zapisz do pliku
    else:
        SeqIO.write( (p1[key])[0], handle_single1, "fastq")
        # Nie ma pary. Zapisz do pliku

for key, value in p2.items():
    if key not in p2.keys():
        SeqIO.write( (p2[key])[0], handle_single2, "fastq")
        # Nie ma pary. Zapisz do pliku
    
handle_paired1.close()
handle_paired2.close()
handle_single1.close()
handle_single2.close()

