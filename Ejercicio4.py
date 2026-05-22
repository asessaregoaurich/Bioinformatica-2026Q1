# coding: utf-8
# !/usr/bin/env python3
"""
Ejercicio 4 - Trabajo Practico Cuatrimestral Parte 2
Bioinformatica - ITBA - 2026Q1

Llama a EMBOSS para:
  1) Calcular ORFs del FASTA de nucleotidos del Ej 1 (getorf) y obtener proteinas
  2) Filtrar solo la proteina del Frame +2 (la confirmada en el Ej 1)
  3) Buscar dominios/motivos PROSITE sobre esa proteina (patmatmotifs)

Output:
  - Ej4_ORFs.fasta            -> todas las proteinas posibles (output de getorf)
  - Ej4_ORF_frame2.fasta      -> solo la proteina del Frame +2
  - Ej4_dominios.patmatmotifs -> reporte crudo de patmatmotifs
  - Ej4_resumen.txt           -> resumen legible de los motivos encontrados
"""

import os
import re
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ============================================================
# CONFIGURACION
# ============================================================

# Input GenBank
INPUT_GENBANK = "sequence.gb"

# FASTA CDS generado automaticamente desde el GB
INPUT_FASTA = "HTT.fasta"

# Outputs
OUT_ORFS_TODOS   = "Ej4_ORFs.fasta"
OUT_DOMINIOS     = "Ej4_dominios.patmatmotifs"
OUT_RESUMEN      = "Ej4_resumen.txt"

# Parametros de getorf
GETORF_FIND    = "1"   # -find 1 -> regiones entre stops, devuelve proteina traducida
GETORF_MINSIZE = "300" # -minsize 300 -> filtra ORFs cortos (HTT esperado ~ 3144 aa)

# Parametros de patmatmotifs
PATMATMOTIFS_FULL = True  # -full -> reporta tambien los motivos no clasicamente fuertes


# ============================================================
# FUNCIONES
# ============================================================

def genbank_a_fasta(genbank_file, output_fasta):
    """
    Extrae TODA la secuencia nucleotidica del GenBank
    y la guarda como FASTA.

    No usa CDS ni coordenadas.
    Exporta el transcript/genoma completo tal como
    aparece en el archivo GenBank.
    """

    if not Path(genbank_file).exists():
        print(f"[ERROR] No se encuentra el archivo GenBank: {genbank_file}")
        sys.exit(1)

    print(f"\n[+] Extrayendo secuencia completa desde {genbank_file}")

    fasta_records = []

    for record in SeqIO.parse(genbank_file, "genbank"):

        gene_name = record.name
        descripcion = record.description

        # TODA la secuencia
        full_seq = record.seq

        fasta_record = SeqRecord(
            full_seq,
            id=gene_name,
            description=descripcion
        )

        fasta_records.append(fasta_record)

    if not fasta_records:
        print("[ERROR] No se pudo leer ninguna secuencia del GenBank.")
        sys.exit(1)

    SeqIO.write(fasta_records, output_fasta, "fasta")

    print(f"[OK] FASTA generado: {output_fasta}")
    print(f"    Secuencias exportadas: {len(fasta_records)}")
    print(f"    Longitud: {len(fasta_records[0].seq)} nt")


def run(cmd, descripcion):
    """Ejecuta un comando externo de EMBOSS y aborta si falla."""
    print(f"\n[+] {descripcion}")
    print(f"    $ {' '.join(cmd)}")
    resultado = subprocess.run(cmd, capture_output=True, text=True)
    if resultado.returncode != 0:
        print(f"[ERROR] Fallo el comando ({resultado.returncode}):")
        print(resultado.stderr)
        sys.exit(1)
    if resultado.stdout.strip():
        print(resultado.stdout)
    return resultado


def verificar_emboss():
    """Chequea que EMBOSS este instalado en el entorno."""
    try:
        subprocess.run(["embossversion", "-stdout", "-auto"],
                       capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[ERROR] EMBOSS no esta instalado o no esta en el PATH.")
        print("        Instalalo con: sudo apt-get install emboss")
        sys.exit(1)


def verificar_prosite():
    """Chequea que PROSITE este disponible para EMBOSS."""
    # 1) Preguntarle a EMBOSS donde estan sus datos
    resultado = subprocess.run(
        ["embossversion", "-full", "-auto"],
        capture_output=True, text=True
    )
    for linea in resultado.stdout.splitlines():
        if "Data" in linea:
            emboss_data_dir = linea.split(":")[1].strip()
            prosite_path = Path(emboss_data_dir) / "PROSITE" / "prosite.lines"
            if prosite_path.exists():
                print(f"[OK] PROSITE encontrado en {prosite_path}")
                return

    # 2) Fallback: variable de entorno EMBOSS_DATA
    emboss_data = os.environ.get("EMBOSS_DATA", "")
    if emboss_data and Path(emboss_data, "prosite.lines").exists():
        print(f"[OK] PROSITE encontrado en {emboss_data}")
        return

    print("[ERROR] No se encuentra el indice PROSITE.")
    print("        Corre: prosextract -prositedir $EMBOSS_DATA")
    sys.exit(1)


def correr_getorf():
    """Llama a getorf para obtener todas las proteinas posibles (los 6 frames)."""
    cmd = [
        "getorf",
        "-sequence", INPUT_FASTA,
        "-outseq",   OUT_ORFS_TODOS,
        "-find",     GETORF_FIND,
        "-minsize",  GETORF_MINSIZE,
        "-auto",
    ]
    run(cmd, f"getorf sobre {INPUT_FASTA} (minsize={GETORF_MINSIZE} aa)")


def correr_patmatmotifs():
    """Corre patmatmotifs sobre cada ORF individualmente."""

    print("\n[+] Ejecutando patmatmotifs ORF por ORF")

    records = list(SeqIO.parse(OUT_ORFS_TODOS, "fasta"))

    salida_total = ""

    for i, rec in enumerate(records, 1):

        temp_fasta = f"tmp_orf_{i}.fasta"

        # guardar ORF individual
        SeqIO.write([rec], temp_fasta, "fasta")

        print(f"    Analizando {rec.id} ({len(rec.seq)} aa)")

        cmd = [
            "patmatmotifs",
            "-sequence", temp_fasta,
            "-stdout",
            "-auto",
        ]

        if PATMATMOTIFS_FULL:
            cmd.insert(1, "-full")

        resultado = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        salida_total += resultado.stdout
        salida_total += "\n\n"

        # borrar temporal
        Path(temp_fasta).unlink()

    # guardar output combinado
    with open(OUT_DOMINIOS, "w") as f:
        f.write(salida_total)

    print(f"[OK] Reporte combinado: {OUT_DOMINIOS}")

def correr_pepstats():
    """Corre pepstats sobre cada ORF individualmente."""

    print("\n[+] Ejecutando pepstats ORF por ORF")

    records = list(SeqIO.parse(OUT_ORFS_TODOS, "fasta"))

    salida_total = ""

    for i, rec in enumerate(records, 1):

        temp_fasta = f"tmp_orf_{i}.fasta"

        # guardar ORF individual
        SeqIO.write([rec], temp_fasta, "fasta")

        print(f"    Analizando {rec.id} ({len(rec.seq)} aa)")

        cmd = [
            "pepstats",
            "-sequence", temp_fasta,
            "-stdout",
            "-auto",
        ]

        resultado = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        salida_total += (
            f"\n{'='*60}\n"
            f"ORF {i}: {rec.id} ({len(rec.seq)} aa)\n"
            f"{'='*60}\n"
        )

        salida_total += resultado.stdout
        salida_total += "\n\n"

        # borrar temporal
        Path(temp_fasta).unlink()

    # guardar output combinado
    with open("Ej4_pepstats.txt", "w") as f:
        f.write(salida_total)

    print("[OK] Reporte combinado: Ej4_pepstats.txt")

def resumir_resultados():
    """Genera un resumen legible de motivos PROSITE por ORF."""

    if not Path(OUT_DOMINIOS).exists():
        print(f"[ERROR] No se encuentra {OUT_DOMINIOS}")
        return

    motivos = []
    secuencia_actual = None
    bloque = {}

    with open(OUT_DOMINIOS, "r") as f:

        for linea in f:
            linea = linea.strip()

            if linea.startswith("# Sequence:"):

                secuencia_actual = (
                    linea.replace("# Sequence:", "")
                    .split("from:")[0]
                    .strip()
                )

            elif linea.startswith("Length"):
                bloque["length"] = linea.split("=")[1].strip()

            elif linea.startswith("Start"):
                bloque["start"] = (
                    linea.split("position")[1]
                    .split("of")[0]
                    .strip()
                )

            elif linea.startswith("End"):
                bloque["end"] = (
                    linea.split("position")[1]
                    .split("of")[0]
                    .strip()
                )

            elif linea.startswith("Motif"):

                bloque["motif"] = linea.split("=")[1].strip()
                bloque["sequence"] = secuencia_actual

                motivos.append(dict(bloque))
                bloque = {}

    # =====================================================
    # ESCRIBIR RESUMEN
    # =====================================================

    with open(OUT_RESUMEN, "w") as f:

        f.write("Ejercicio 4 - Resumen de motivos PROSITE\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Input proteinas : {OUT_ORFS_TODOS}\n")
        f.write(f"Reporte crudo   : {OUT_DOMINIOS}\n")
        f.write(f"Motivos hallados: {len(motivos)}\n\n")

        por_orf = {}

        for m in motivos:
            seq = m["sequence"]

            if seq not in por_orf:
                por_orf[seq] = []

            por_orf[seq].append(m)

        for orf, lista in por_orf.items():

            f.write("=" * 60 + "\n")
            f.write(f"ORF: {orf}\n")
            f.write(f"Motivos encontrados: {len(lista)}\n")
            f.write("=" * 60 + "\n\n")

            for i, m in enumerate(lista, 1):

                f.write(f"[{i}] Motif: {m['motif']}\n")
                f.write(f"    Posicion: {m['start']} - {m['end']}\n")
                f.write(f"    Longitud: {m['length']}\n\n")

    print(f"\n[OK] Resumen escrito en {OUT_RESUMEN}")

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("Ejercicio 4 - EMBOSS + PROSITE")
    print("=" * 60)

    if not Path(INPUT_GENBANK).exists():
        print(f"[ERROR] No se encuentra el GenBank de entrada: {INPUT_GENBANK}")
        sys.exit(1)

    verificar_emboss()
    verificar_prosite()

    #GB -> FASTA
    genbank_a_fasta(INPUT_GENBANK, INPUT_FASTA)

    #Pipeline principal
    correr_getorf()
    correr_patmatmotifs()
    correr_pepstats()
    resumir_resultados()

    print("\n[OK] Ejercicio 4 finalizado.")


if __name__ == "__main__":
    main()