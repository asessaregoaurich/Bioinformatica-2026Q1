# coding: utf-8
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# 1. Lectura del archivo GenBank 
input_file = "sequence.gb"
record = SeqIO.read(input_file, "genbank")
dna_seq = record.seq

# 2. Generación de los 6 marcos de lectura (Reading Frames) 
# 3 marcos de la cadena directa y 3 de la complementaria reversa
frames = []

# Cadena Directa (Sense)
for i in range(3):
    frames.append({
        "id": f"Frame_+{i+1}",
        "seq": dna_seq[i:]
    })

# Cadena Complementaria Reversa (Antisense)
rev_seq = dna_seq.reverse_complement()
for i in range(3):
    frames.append({
        "id": f"Frame_-{i+1}",
        "seq": rev_seq[i:]
    })

# 3. Traducción y creación de registros FASTA [cite: 34, 39]
protein_records = []

for frame in frames:
    # Traducimos la secuencia de nucleótidos a aminoácidos
    # El parámetro to_stop=False permite ver toda la secuencia (incluyendo '*' para stops)
    protein_seq = frame["seq"].translate()
    
    # Creamos un objeto SeqRecord para cada traducción
    new_record = SeqRecord(
        protein_seq, 
        id=f"{record.id}_{frame['id']}", 
        description=f"Traduccion automatica del gen {record.name}"
    )
    protein_records.append(new_record)

# 4. Escritura del archivo de salida en formato FASTA [cite: 34, 39]
output_file = "Ej1_ORF.fas"
SeqIO.write(protein_records, output_file, "fasta")

print(f"Proceso completado. Se han generado 6 traducciones en '{output_file}'.")