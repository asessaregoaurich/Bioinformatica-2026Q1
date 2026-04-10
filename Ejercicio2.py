#!/usr/bin/env python3
"""
Ejercicio 2.a - BLAST local y/o remoto
Modificá las variables de configuración abajo antes de correr.
"""

import os
import sys
import time
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# ============================================================
#  CONFIGURACIÓN — modificá estas rutas
# ============================================================

INPUT_FASTA = "Ej1_ORF.fas"            # archivo FASTA de aminoácidos
DB_LOCAL    = r"C:\Users\acorb\Documents\ITBA\1C_2026\Bioinformática\Bioinformatica-2026Q1\swissprot"    # ruta a la DB local (sin extensión)
OUTPUT_DIR  = "blast_results"          # carpeta donde se guardan los resultados

CORRER_REMOTO = False   # True para correr BLAST remoto
CORRER_LOCAL  = True    # True para correr BLAST local

# ============================================================


def blast_remoto(seq_record, output_dir):
    seq_id = seq_record.id.replace("|", "_").replace("/", "_")
    output_file = os.path.join(output_dir, f"blast_remote_{seq_id}.xml")

    print(f"  [REMOTO] Ejecutando BLASTp para: {seq_record.id}")
    print(f"  Largo: {len(seq_record.seq)} aa — enviando a NCBI (puede tardar)...")

    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="swissprot",
        sequence=str(seq_record.seq),
        hitlist_size=10,
        expect=0.001,
        format_type="XML"
    )

    with open(output_file, "w") as out:
        out.write(result_handle.read())
    result_handle.close()

    print(f"  Guardado en: {output_file}")
    return output_file


def blast_local(seq_record, db_path, output_dir):
    seq_id = seq_record.id.replace("|", "_").replace("/", "_")
    output_file = os.path.join(output_dir, f"blast_local_{seq_id}.xml")

    print(f"  [LOCAL] Ejecutando BLASTp para: {seq_record.id}")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">{seq_record.id}\n{str(seq_record.seq)}\n")
        tmp_path = tmp.name

    try:
        cmd = [
            "blastp",
            "-query", tmp_path,
            "-db", db_path,
            "-out", output_file,
            "-outfmt", "5",
            "-max_target_seqs", "10",
            "-evalue", "0.001"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  ERROR: {result.stderr}")
            return None

        print(f"  Guardado en: {output_file}")
        return output_file
    finally:
        os.unlink(tmp_path)


def mostrar_resumen(xml_file, modo):
    print(f"\n  --- Resumen {modo} ---")
    with open(xml_file) as f:
        blast_records = list(NCBIXML.parse(f))

    tiene_hits = False
    for blast_record in blast_records:
        print(f"  Query: {blast_record.query}")
        print(f"  Hits encontrados: {len(blast_record.alignments)}\n")

        if not blast_record.alignments:
            print("  No se encontraron hits significativos.")
            continue

        tiene_hits = True
        print(f"  {'#':<4} {'Descripción':<50} {'Score':>8} {'E-value':>12} {'Identity':>10}")
        print(f"  {'-'*88}")

        for i, alignment in enumerate(blast_record.alignments):
            hsp = alignment.hsps[0]
            identity_pct = (hsp.identities / hsp.align_length) * 100
            desc = alignment.title[:48] + ".." if len(alignment.title) > 48 else alignment.title
            print(f"  {i+1:<4} {desc:<50} {hsp.score:>8.0f} {hsp.expect:>12.2e} {identity_pct:>9.1f}%")
        print()

    # Borrar el archivo si no tuvo hits
    if not tiene_hits:
        os.remove(xml_file)
        print(f"  Archivo eliminado (sin hits): {xml_file}")


def main():
    if not os.path.exists(INPUT_FASTA):
        print(f"ERROR: No se encontró el archivo: {INPUT_FASTA}")
        sys.exit(1)

    if CORRER_LOCAL and not os.path.exists(DB_LOCAL + ".pin") and not os.path.exists(DB_LOCAL + ".psq"):
        print(f"ERROR: No se encontró la DB local en '{DB_LOCAL}'.")
        print("Asegurate de que la carpeta swissprot esté en el mismo lugar que este script.")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    sequences = list(SeqIO.parse(INPUT_FASTA, "fasta"))
    if not sequences:
        print(f"ERROR: No se encontraron secuencias en {INPUT_FASTA}")
        sys.exit(1)

    print(f"\n=== Ex2 - BLAST ===")
    print(f"Input: {INPUT_FASTA} ({len(sequences)} secuencias)")
    print(f"Output: {OUTPUT_DIR}/\n")

    for i, seq_record in enumerate(sequences):
        print(f"\n[{i+1}/{len(sequences)}] {seq_record.id}")

        if CORRER_REMOTO:
            xml = blast_remoto(seq_record, OUTPUT_DIR)
            if xml:
                mostrar_resumen(xml, "REMOTO")
            if i < len(sequences) - 1:
                print("  Esperando 15s (límite NCBI)...")
                time.sleep(15)

        if CORRER_LOCAL:
            xml = blast_local(seq_record, DB_LOCAL, OUTPUT_DIR)
            if xml:
                mostrar_resumen(xml, "LOCAL")

    print(f"\n=== Listo. Resultados en: {OUTPUT_DIR}/ ===")


if __name__ == "__main__":
    main()

