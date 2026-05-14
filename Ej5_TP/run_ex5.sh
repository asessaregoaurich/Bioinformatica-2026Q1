#!/bin/bash
# ============================================================
#  run_ex5.sh  –  Ejercicio 5: Diseño de Primers
#  Introducción a la Bioinformática - TP Cuatrimestral
#
#  Uso:
#      bash run_ex5.sh <transcripto.fasta> [primer_config.json]
#
#  Ejemplo:
#      bash run_ex5.sh HTT_transcripto.fasta primer_config.json
# ============================================================

set -euo pipefail

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

INPUT_FASTA="${1:-}"
CONFIG="${2:-primer_config.json}"
OUTPUT_DIR="ex5_output"
OUTPUT_FILE="${OUTPUT_DIR}/primers_output.txt"
LOG_FILE="${OUTPUT_DIR}/ex5_run.log"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

log() { echo -e "$1" | tee -a "$LOG_FILE"; }
error_exit() { log "${RED}[ERROR] $1${NC}"; exit 1; }

# ── Validaciones ────────────────────────────────────────────
[[ -z "$INPUT_FASTA" ]] && { echo "Uso: bash run_ex5.sh <fasta> [config.json]"; exit 1; }
[[ ! -f "$INPUT_FASTA" ]] && error_exit "No se encontró: $INPUT_FASTA"
[[ ! -f "$CONFIG" ]]      && error_exit "No se encontró config: $CONFIG"

mkdir -p "$OUTPUT_DIR"

log ""
log "${BLUE}============================================================${NC}"
log "${BLUE}  EJERCICIO 5 - Diseño de Primers${NC}"
log "${BLUE}  Iniciado: ${TIMESTAMP}${NC}"
log "${BLUE}============================================================${NC}"
log ""
log "FASTA entrada : $INPUT_FASTA"
log "Config JSON   : $CONFIG"
log "Salida        : $OUTPUT_FILE"
log ""

# ── BioPython ───────────────────────────────────────────────
if ! python3 -c "import Bio" &>/dev/null; then
    log "${YELLOW}[AVISO]${NC} Instalando BioPython..."
    pip3 install biopython --quiet
fi

# ── Ejecutar ─────────────────────────────────────────────────
python3 Ex5.py -i "$INPUT_FASTA" -c "$CONFIG" -o "$OUTPUT_FILE" 2>&1 | tee -a "$LOG_FILE"

log ""
log "${GREEN}[COMPLETADO]${NC} Primers guardados en: $OUTPUT_FILE"
log "Log: $LOG_FILE"
