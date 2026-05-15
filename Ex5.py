#!/usr/bin/env python3
"""
Ejercicio 5 - Diseño de Primers
Introducción a la Bioinformática - TP Cuatrimestral
Gen: HTT (Huntingtina) - Enfermedad de Huntington

Descripción:
    Lee una secuencia de transcripto en formato FASTA y diseña primers
    que cumplan los criterios definidos en un archivo de configuración JSON.

    Criterios (configurables en primer_config.json):
        - Longitud: 18 a 24 pares de bases
        - Contenido GC: 50% a 60%
        - Sin GC en el extremo 3'
        - Temperatura de melting <= 67°C

Uso:
    python Ex5.py -i transcripto.fasta -c primer_config.json [-o primers_output.txt]

Requisitos:
    - BioPython: pip install biopython
"""

import argparse
import json
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq


# ─────────────────────────────────────────────────────────────
# Cálculo de Temperatura de Melting
# ─────────────────────────────────────────────────────────────

def calc_tm(sequence: str) -> float:
    """
    Calcula la temperatura de melting (Tm) de un primer.

    Usa dos fórmulas según el largo del primer:
      - Para primers cortos (<= 13 nt): Wallace Rule
            Tm = 2*(A+T) + 4*(G+C)
      - Para primers más largos (> 13 nt): fórmula de nearest-neighbor simplificada
            Tm = 64.9 + 41*(G+C - 16.4) / longitud

    Ambas son aproximaciones estándar en diseño de primers.
    """
    seq = sequence.upper()
    length = len(seq)
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')

    if length <= 13:
        # Wallace Rule
        tm = 2 * at + 4 * gc
    else:
        # Fórmula simplificada de nearest-neighbor
        tm = 64.9 + 41 * (gc - 16.4) / length

    return round(tm, 2)


# ─────────────────────────────────────────────────────────────
# Cálculo de contenido GC
# ─────────────────────────────────────────────────────────────

def calc_gc(sequence: str) -> float:
    """Retorna el porcentaje de GC de una secuencia."""
    seq = sequence.upper()
    gc = seq.count('G') + seq.count('C')
    return round((gc / len(seq)) * 100, 2)


# ─────────────────────────────────────────────────────────────
# Validación de un candidato a primer
# ─────────────────────────────────────────────────────────────

def is_valid_primer(sequence: str, config: dict) -> tuple[bool, dict]:
    """
    Evalúa si una secuencia cumple todos los criterios del config.

    Retorna:
        (True, métricas) si cumple todos los criterios
        (False, métricas) si falla alguno, con métricas igual calculadas
    """
    cfg = config["primer_design"]
    seq = sequence.upper()
    length = len(seq)

    gc_percent = calc_gc(seq)
    tm = calc_tm(seq)
    ends_in_gc = seq[-1] in ('G', 'C')
    starts_in_gc = seq[0] in ('G', 'C')

    metrics = {
        "sequence":   seq,
        "length":     length,
        "gc_percent": gc_percent,
        "tm":         tm,
        "ends_in_gc": ends_in_gc,
    }

    # ── Aplicar criterios ──────────────────────────────────
    if not (cfg["length"]["min"] <= length <= cfg["length"]["max"]):
        return False, metrics

    if not (cfg["gc_content"]["min_percent"] <= gc_percent <= cfg["gc_content"]["max_percent"]):
        return False, metrics

    if tm > cfg["melting_temperature"]["max_celsius"]:
        return False, metrics

    if cfg["gc_clamp"]["avoid_gc_at_3prime"] and ends_in_gc:
        return False, metrics

    if cfg["gc_clamp"]["avoid_gc_at_5prime"] and starts_in_gc:
        return False, metrics

    return True, metrics


# ─────────────────────────────────────────────────────────────
# Diseño de primers por ventana deslizante
# ─────────────────────────────────────────────────────────────

def design_primers(sequence: str, config: dict) -> list[dict]:
    """
    Recorre la secuencia con una ventana deslizante probando todos los
    posibles primers de largo entre min y max definidos en el config.

    Para cada posición y largo, evalúa tanto la cadena directa (forward)
    como el complemento reverso (reverse).

    Retorna una lista de primers válidos ordenados por Tm descendente.
    """
    cfg = config["primer_design"]
    seq = sequence.upper()
    seq_len = len(seq)
    step = cfg.get("slide_step", 1)

    candidates = []
    seen = set()  # evitar duplicados

    min_len = cfg["length"]["min"]
    max_len = cfg["length"]["max"]

    for start in range(0, seq_len - min_len, step):
        for plen in range(min_len, max_len + 1):
            end = start + plen
            if end > seq_len:
                break

            # ── Forward primer ──────────────────────────
            fwd_seq = seq[start:end]
            if fwd_seq not in seen:
                valid, metrics = is_valid_primer(fwd_seq, config)
                if valid:
                    candidates.append({
                        **metrics,
                        "direction": "forward",
                        "start":     start + 1,   # posición 1-based
                        "end":       end,
                    })
                    seen.add(fwd_seq)

            # ── Reverse primer (complemento reverso) ────
            rev_seq = str(Seq(fwd_seq).reverse_complement())
            if rev_seq not in seen:
                valid, metrics = is_valid_primer(rev_seq, config)
                if valid:
                    candidates.append({
                        **metrics,
                        "direction": "reverse",
                        "start":     start + 1,
                        "end":       end,
                    })
                    seen.add(rev_seq)

    # Ordenar por Tm descendente (primers con Tm más alta primero)
    candidates.sort(key=lambda x: x["tm"], reverse=True)

    return candidates


# ─────────────────────────────────────────────────────────────
# Selección de los N mejores primers (diversidad de posiciones)
# ─────────────────────────────────────────────────────────────

def select_best_primers(candidates: list[dict], n: int, min_distance: int = 50) -> list[dict]:
    """
    Selecciona los N mejores primers intentando que estén distribuidos
    a lo largo de la secuencia (separados al menos min_distance nt entre sí).

    Estrategia: greedy — toma el mejor candidato disponible, luego descarta
    los que estén demasiado cerca, y repite.
    """
    selected = []
    used_positions = []

    for candidate in candidates:
        if len(selected) >= n:
            break

        pos = candidate["start"]

        # Verificar que no esté demasiado cerca de los ya seleccionados
        too_close = any(abs(pos - p) < min_distance for p in used_positions)
        if not too_close:
            selected.append(candidate)
            used_positions.append(pos)

    # Si no hubo suficientes con min_distance, completar sin esa restricción
    if len(selected) < n:
        for candidate in candidates:
            if len(selected) >= n:
                break
            if candidate not in selected:
                selected.append(candidate)

    return selected


# ─────────────────────────────────────────────────────────────
# Escritura de resultados
# ─────────────────────────────────────────────────────────────

def write_output(primers: list[dict], output_file: str, config: dict, gene_name: str):
    """Escribe los primers seleccionados en un archivo de texto formateado."""

    cfg = config["primer_design"]

    with open(output_file, "w") as f:
        f.write("=" * 65 + "\n")
        f.write("  DISEÑO DE PRIMERS - EJERCICIO 5\n")
        f.write(f"  Gen: {gene_name}\n")
        f.write("=" * 65 + "\n\n")

        f.write("PARÁMETROS DE DISEÑO (desde primer_config.json):\n")
        f.write(f"  Longitud          : {cfg['length']['min']} - {cfg['length']['max']} nt\n")
        f.write(f"  Contenido GC      : {cfg['gc_content']['min_percent']}% - {cfg['gc_content']['max_percent']}%\n")
        f.write(f"  Tm máxima         : {cfg['melting_temperature']['max_celsius']}°C\n")
        f.write(f"  Evitar GC en 3'   : {cfg['gc_clamp']['avoid_gc_at_3prime']}\n")
        f.write(f"  Cantidad pedida   : {cfg['num_primers']} primers\n")
        f.write("\n" + "-" * 65 + "\n\n")

        if not primers:
            f.write("No se encontraron primers que cumplan todos los criterios.\n")
            f.write("Sugerencia: relajar los parámetros en primer_config.json\n")
            return

        for i, p in enumerate(primers, 1):
            f.write(f"PRIMER #{i}\n")
            f.write(f"  Secuencia   : 5'-{p['sequence']}-3'\n")
            f.write(f"  Dirección   : {p['direction']}\n")
            f.write(f"  Posición    : {p['start']} - {p['end']} (en el transcripto)\n")
            f.write(f"  Longitud    : {p['length']} nt\n")
            f.write(f"  %GC         : {p['gc_percent']}%\n")
            f.write(f"  Tm          : {p['tm']}°C\n")
            f.write(f"  Termina en GC: {'Sí' if p['ends_in_gc'] else 'No'}\n")
            f.write("\n")

        f.write("-" * 65 + "\n")
        f.write(f"Total de primers diseñados: {len(primers)}\n")


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Ejercicio 5 - Diseño de primers para gen HTT (Huntington)"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Archivo FASTA con la secuencia del transcripto"
    )
    parser.add_argument(
        "-c", "--config",
        default="primer_config.json",
        help="Archivo de configuración JSON (default: primer_config.json)"
    )
    parser.add_argument(
        "-o", "--output",
        default="primers_output.txt",
        help="Archivo de salida con los primers (default: primers_output.txt)"
    )

    args = parser.parse_args()

    # ── Validaciones ────────────────────────────────────────
    if not os.path.exists(args.input):
        print(f"[ERROR] No se encontró el archivo: {args.input}")
        sys.exit(1)

    if not os.path.exists(args.config):
        print(f"[ERROR] No se encontró el config: {args.config}")
        sys.exit(1)

    # ── Leer configuración ──────────────────────────────────
    with open(args.config, "r") as f:
        config = json.load(f)

    gene_name = config.get("gene_info", {}).get("name", "desconocido")
    n_primers = config["primer_design"]["num_primers"]

    print("=" * 65)
    print(f"  EJERCICIO 5 - Diseño de Primers")
    print(f"  Gen: {gene_name}  |  Config: {args.config}")
    print("=" * 65)

    # ── Leer secuencia ──────────────────────────────────────
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        print("[ERROR] No se encontraron secuencias en el archivo FASTA.")
        sys.exit(1)

    # Usar la primera secuencia del archivo
    record = records[0]
    sequence = str(record.seq).upper().replace("\n", "")
    print(f"\nSecuencia cargada : {record.id}")
    print(f"Longitud          : {len(sequence)} nt")

    # ── Diseñar primers ─────────────────────────────────────
    print(f"\nBuscando candidatos en los {len(sequence)} nucleótidos...")
    all_candidates = design_primers(sequence, config)
    print(f"Candidatos válidos encontrados: {len(all_candidates)}")

    if not all_candidates:
        print("\n[AVISO] No se encontraron primers con los criterios actuales.")
        print("  Sugerencias:")
        print("  - Aumentar rango de longitud en el JSON")
        print("  - Ampliar rango de %GC")
        print("  - Aumentar Tm máxima")
        write_output([], args.output, config, gene_name)
        sys.exit(0)

    # ── Seleccionar los mejores ─────────────────────────────
    best = select_best_primers(all_candidates, n_primers)
    print(f"Primers seleccionados: {len(best)}")

    # ── Mostrar resumen en consola ──────────────────────────
    print("\n" + "-" * 65)
    print(f"{'#':<4} {'Secuencia':<26} {'Dir':<8} {'Long':>4} {'GC%':>6} {'Tm':>7}  Pos")
    print("-" * 65)
    for i, p in enumerate(best, 1):
        print(
            f"{i:<4} {p['sequence']:<26} {p['direction']:<8} "
            f"{p['length']:>4} {p['gc_percent']:>5}% {p['tm']:>6}°C  "
            f"{p['start']}-{p['end']}"
        )
    print("-" * 65)

    # ── Escribir archivo de salida ──────────────────────────
    write_output(best, args.output, config, gene_name)
    print(f"\n[OK] Resultados guardados en: {args.output}")


if __name__ == "__main__":
    main()
