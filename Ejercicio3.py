# coding: utf-8
from Bio import Entrez, SeqIO, SearchIO
import os

# 1. Configuración de seguridad y contacto para NCBI
# NCBI pide un mail para contactarte si tus consultas saturan sus servidores
Entrez.email = "jaldana@itba.edu.ar" 

# 2. Definición de archivos de entrada y salida
blast_xml = "blast_results/blast_local_NM_001388492.1_Frame_+2.xml" # El output del Ejercicio 2
output_fasta = "msa_input.fas"          # El archivo que contendrá todas las secuencias
query_fasta = "Ej1_ORF.fas"            # Tu secuencia original (Frame +2)

def generar_msa_input():
    print("Iniciando la preparación de datos para el MSA...")
    
    # Lista para guardar los registros de proteínas que vamos a alinear
    sequences_to_align = []

    # 3. Cargar nuestra secuencia original (Query)
    # Solo tomamos el Frame +2, que es el que validamos como correcto
    for record in SeqIO.parse(query_fasta, "fasta"):
        if "Frame_+2" in record.id:
            record.id = "Human_HTT_Query"
            record.description = "Secuencia consulta - Homo sapiens"
            sequences_to_align.append(record)
            break

    # 4. Parsear el resultado del BLAST para obtener los 10 mejores hits
    if not os.path.exists(blast_xml):
        print(f"Error: No se encuentra el archivo {blast_xml}. ¿Corriste el Ejercicio 2?")
        return

    # Usamos SearchIO para leer el XML de BLAST de forma sencilla
    blast_qresult = SearchIO.read(blast_xml, 'blast-xml')
    
    # Tomamos los IDs de los primeros 10 hits (excluyendo el hit 1 si es idéntico a la query)
    hit_ids = []
    for hit in blast_qresult[:11]: # Pedimos 11 por si el primero es nuestra propia secuencia
        # El ID suele venir como 'sp|P42858|HD_HUMAN', extraemos solo el ID de acceso 'P42858'
        accession = hit.id.split('|')[1] if '|' in hit.id else hit.id
        if accession not in [s.id for s in sequences_to_align]:
            hit_ids.append(accession)
    
    # Limitamos a 10 hits según pide el enunciado
    hit_ids = hit_ids[:10]

    # 5. Descargar las secuencias desde la base de datos de Proteínas del NCBI
    print(f"Descargando {len(hit_ids)} secuencias desde NCBI/UniProt...")
    for uid in hit_ids:
        try:
            # efetch descarga el registro completo en formato FASTA
            handle = Entrez.efetch(db="protein", id=uid, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            sequences_to_align.append(record)
            handle.close()
            print(f"  OK: Descargado {uid}")
        except Exception as e:
            print(f"  Error al descargar {uid}: {e}")

    # 6. Guardar todo en un único archivo FASTA
    SeqIO.write(sequences_to_align, output_fasta, "fasta")
    print(f"\nListo! El archivo '{output_fasta}' fue generado con {len(sequences_to_align)} secuencias.")

if __name__ == "__main__":
    generar_msa_input()