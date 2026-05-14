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

# ============================================================
# CONFIGURACION
# ============================================================

# Input: FASTA de nucleotidos del Ej 1
INPUT_FASTA = "HTT_NM_001388492.fasta"

# Outputs
OUT_ORFS_TODOS   = "Ej4_ORFs.fasta"
OUT_ORF_FRAME2   = "Ej4_ORF_frame2.fasta"
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


def extraer_frame2():
    """
    getorf rotula cada ORF con la posicion de inicio y la hebra.
    Para sentido directo: [start - end]   -> Frame +1, +2 o +3 segun start mod 3
    Para hebra reversa:   [end - start] (REVERSE SENSE)

    Frame +2 == start (1-indexado) tiene resto 2 al dividir por 3.
    Tomamos el ORF mas largo del Frame +2 (asumimos que ese es el del HTT).
    """
    if not Path(OUT_ORFS_TODOS).exists():
        print(f"[ERROR] No se encuentra {OUT_ORFS_TODOS}")
        sys.exit(1)

    print(f"\n[+] Filtrando ORFs del Frame +2 desde {OUT_ORFS_TODOS}")

    # Parsear el FASTA con BioPython
    records = list(SeqIO.parse(OUT_ORFS_TODOS, "fasta"))
    print(f"    Total de ORFs leidos: {len(records)}")

    # Buscar candidatos del Frame +2 (hebra directa, start mod 3 == 2)
    candidatos = []
    patron = re.compile(r"\[(\d+)\s*-\s*(\d+)\]")

    for rec in records:
        # rec.description contiene el header completo (incluye [start - end])
        if "REVERSE SENSE" in rec.description:
            continue
        m = patron.search(rec.description)
        if not m:
            continue
        start = int(m.group(1))
        # Frame +1 -> start%3 == 1 ; Frame +2 -> start%3 == 2 ; Frame +3 -> start%3 == 0
        if start % 3 == 2:
            candidatos.append(rec)

    if not candidatos:
        print("[ERROR] No se encontraron ORFs en Frame +2. Revisa el FASTA o bajar minsize.")
        sys.exit(1)

    # Tomar el ORF mas largo (es el ORF principal del HTT)
    rec_elegido = max(candidatos, key=lambda r: len(r.seq))
    print(f"    -> Frame +2 elegido: {rec_elegido.description}  ({len(rec_elegido.seq)} aa)")

    # Escribir como FASTA estandar (BioPython parte en lineas de 60 automaticamente)
    SeqIO.write([rec_elegido], OUT_ORF_FRAME2, "fasta")


def correr_patmatmotifs():
    """Llama a patmatmotifs sobre el ORF del Frame +2."""
    cmd = [
        "patmatmotifs",
        "-sequence", OUT_ORF_FRAME2,
        "-outfile",  OUT_DOMINIOS,
        "-auto",
    ]
    if PATMATMOTIFS_FULL:
        cmd.insert(1, "-full")
    run(cmd, f"patmatmotifs sobre {OUT_ORF_FRAME2}")


def resumir_resultados():
    """Genera Ej4_resumen.txt listando los motivos PROSITE encontrados."""
    if not Path(OUT_DOMINIOS).exists():
        print(f"[ERROR] No se encuentra {OUT_DOMINIOS}")
        return

    motivos = []
    bloque  = {}
    with open(OUT_DOMINIOS) as f:
        for linea in f:
            linea = linea.rstrip()
            if linea.startswith("Length"):
                bloque["length"] = linea.split("=")[1].strip()
            elif linea.startswith("Start"):
                bloque["start"] = linea.split("=")[1].strip()
            elif linea.startswith("End"):
                bloque["end"] = linea.split("=")[1].strip()
            elif linea.startswith("Motif"):
                bloque["motif"] = linea.split("=")[1].strip()
                motivos.append(dict(bloque))
                bloque = {}

    with open(OUT_RESUMEN, "w") as f:
        f.write("Ejercicio 4 - Resumen de motivos PROSITE encontrados\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Input proteina  : {OUT_ORF_FRAME2}\n")
        f.write(f"Reporte crudo   : {OUT_DOMINIOS}\n")
        f.write(f"Motivos hallados: {len(motivos)}\n\n")
        for i, m in enumerate(motivos, 1):
            f.write(f"[{i}] Motif: {m.get('motif','?')}\n")
            f.write(f"    Posicion: {m.get('start','?')} - {m.get('end','?')}\n")
            f.write(f"    Longitud del match: {m.get('length','?')}\n\n")

    print(f"\n[OK] Resumen escrito en {OUT_RESUMEN} ({len(motivos)} motivos)")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("Ejercicio 4 - EMBOSS + PROSITE")
    print("=" * 60)

    if not Path(INPUT_FASTA).exists():
        print(f"[ERROR] No se encuentra el FASTA de entrada: {INPUT_FASTA}")
        sys.exit(1)

    verificar_emboss()
    verificar_prosite()
    correr_getorf()
    extraer_frame2()
    correr_patmatmotifs()
    resumir_resultados()

    print("\n[OK] Ejercicio 4 finalizado.")


if __name__ == "__main__":
    main()