#!/bin/bash
# =============================================================
# run.sh — Pipeline TP Cuatrimestral Bioinformática ITBA 2026Q1
# =============================================================

source ~/bioenv/bin/activate

# =============================================================
# LOGGING
# =============================================================

LOG_FILE="pipeline.log"

log() {
    local nivel="$1"
    local mensaje="$2"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${nivel} - ${mensaje}" | tee -a "$LOG_FILE"
}

separador() {
    echo "================================================================" | tee -a "$LOG_FILE"
}

# =============================================================
# USO
# =============================================================

usage() {
    echo "Uso: $0 {1|2|3|4|5|all} [args...]"
    echo "  1:   Procesamiento de Secuencias"
    echo "  2:   BLAST"
    echo "  3:   Alineamiento Múltiple - MSA"
    echo "  4:   Dominios y Motivos PROSITE"
    echo "  5:   Diseño de Primers [fasta] [config.json]"
    echo "  all: Corre toda la secuencia (ej5 usa inputs por defecto)"
    exit 1
}

if [ -z "$1" ]; then usage; fi

# =============================================================
# EJERCICIO 1
# =============================================================

ejercicio_1() {
    separador
    log "INFO" "Iniciando Ejercicio 1 — Procesamiento de Secuencias"

    # Validación de input
    if [ ! -f "inputs/sequence.gb" ]; then
        log "ERROR" "No se encontró inputs/sequence.gb"
        exit 1
    fi
    log "INFO" "Input validado: inputs/sequence.gb"

    # Ejecución
    python3 TP_Parte1/Ej1/ejercicio1.py >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR" "Ejercicio 1 falló. Revisar log para detalles."
        exit 1
    fi

    # Validación de output
    if [ ! -f "TP_Parte1/Ej1/output/Ej1_ORF.fas" ]; then
        log "ERROR" "Output no generado: TP_Parte1/Ej1/output/Ej1_ORF.fas"
        exit 1
    fi
    log "OK" "Ejercicio 1 finalizado. Output: TP_Parte1/Ej1/output/Ej1_ORF.fas"
}

# =============================================================
# EJERCICIO 2
# =============================================================

ejercicio_2() {
    separador
    log "INFO" "Iniciando Ejercicio 2 — BLAST"

    # Validación de input
    if [ ! -f "TP_Parte1/Ej1/output/Ej1_ORF.fas" ]; then
        log "ERROR" "No se encontró TP_Parte1/Ej1/output/Ej1_ORF.fas. ¿Corriste el Ejercicio 1?"
        exit 1
    fi
    log "INFO" "Input validado: TP_Parte1/Ej1/output/Ej1_ORF.fas"

    # Ejecución
    python3 TP_Parte1/Ej2/Ejercicio2.py >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR" "Ejercicio 2 falló. Revisar log para detalles."
        exit 1
    fi

    # Validación de output (al menos un XML generado)
    xml_count=$(ls TP_Parte1/Ej2/output/*.xml 2>/dev/null | wc -l)
    if [ "$xml_count" -eq 0 ]; then
        log "ERROR" "No se generaron archivos XML en TP_Parte1/Ej2/output/"
        exit 1
    fi
    log "OK" "Ejercicio 2 finalizado. Archivos XML generados: $xml_count"
}

# =============================================================
# EJERCICIO 3
# =============================================================

ejercicio_3() {
    separador
    log "INFO" "Iniciando Ejercicio 3 — MSA"

    # Validación de inputs
    if [ ! -f "TP_Parte1/Ej2/output/blast_local_NM_001388492.1_Frame_+2.xml" ]; then
        log "ERROR" "No se encontró blast_local_NM_001388492.1_Frame_+2.xml. ¿Corriste el Ejercicio 2 en modo local?"
        exit 1
    fi
    if [ ! -f "TP_Parte1/Ej1/output/Ej1_ORF.fas" ]; then
        log "ERROR" "No se encontró Ej1_ORF.fas. ¿Corriste el Ejercicio 1?"
        exit 1
    fi
    log "INFO" "Inputs validados: blast XML + Ej1_ORF.fas"

    # Validación de MUSCLE
    if ! command -v muscle &> /dev/null; then
    log "ERROR" "MUSCLE no está instalado. Instalá con: sudo apt install muscle"
    exit 1
    fi
    
    # Ejecución
    python3 TP_Parte1/Ej3/Ejercicio3.py >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR" "Ejercicio 3 falló. Revisar log para detalles."
        exit 1
    fi

    # Validación de output
    if [ ! -f "TP_Parte1/Ej3/output/msa_output.aln" ]; then
        log "ERROR" "Output no generado: TP_Parte1/Ej3/output/msa_output.aln"
        exit 1
    fi
    log "OK" "Ejercicio 3 finalizado. Output: TP_Parte1/Ej3/output/msa_output.aln"
}

# =============================================================
# EJERCICIO 4
# =============================================================

ejercicio_4() {
    separador
    log "INFO" "Iniciando Ejercicio 4 — EMBOSS + PROSITE"

    # Validación de inputs
    if [ ! -f "inputs/sequence.gb" ]; then
        log "ERROR" "No se encontró inputs/sequence.gb"
        exit 1
    fi
    log "INFO" "Input validado: inputs/sequence.gb"

    # Validación de EMBOSS
    if ! command -v patmatmotifs &> /dev/null; then
        log "ERROR" "EMBOSS no está instalado. Instalá con: sudo apt install emboss"
        exit 1
    fi
    log "INFO" "EMBOSS verificado"

    # Indexado de PROSITE si es necesario
    if [ ! -f "prosite_data/PROSITE/prosite.lines" ]; then
        log "INFO" "Indexando PROSITE por primera vez..."
        mkdir -p prosite_data/PROSITE

        if [ ! -f "prosite_data/PROSITE/prosite.dat" ]; then
            log "INFO" "Descargando prosite.dat..."
            wget https://ftp.expasy.org/databases/prosite/prosite.dat -P prosite_data/PROSITE/ >> "$LOG_FILE" 2>&1
            if [ $? -ne 0 ]; then
                log "ERROR" "Falló la descarga de prosite.dat"
                exit 1
            fi
        fi

        if [ ! -f "prosite_data/PROSITE/prosite.doc" ]; then
            log "INFO" "Descargando prosite.doc..."
            wget https://ftp.expasy.org/databases/prosite/prosite.doc -P prosite_data/PROSITE/ >> "$LOG_FILE" 2>&1
            if [ $? -ne 0 ]; then
                log "ERROR" "Falló la descarga de prosite.doc"
                exit 1
            fi
        fi

        export EMBOSS_DATA=prosite_data/
        prosextract -prositedir prosite_data/PROSITE/ >> "$LOG_FILE" 2>&1
        if [ $? -ne 0 ]; then
            log "ERROR" "Falló prosextract"
            exit 1
        fi
        log "INFO" "PROSITE indexado correctamente"
    else
        export EMBOSS_DATA=prosite_data/
        log "INFO" "PROSITE ya indexado, saltando prosextract"
    fi

    # Ejecución
    python3 TP_Parte2/Ej4/Ejercicio4.py >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR" "Ejercicio 4 falló. Revisar log para detalles."
        exit 1
    fi

    # Validación de outputs
    if [ ! -f "TP_Parte2/Ej4/output/Ej4_resumen.txt" ]; then
        log "ERROR" "Output no generado: TP_Parte2/Ej4/output/Ej4_resumen.txt"
        exit 1
    fi
    log "OK" "Ejercicio 4 finalizado. Outputs en: TP_Parte2/Ej4/output/"
}

# =============================================================
# EJERCICIO 5
# =============================================================

ejercicio_5() {
    separador
    local input_fasta="${1:-inputs/HTT.fasta}"
    local config="${2:-inputs/primer_config.json}"

    log "INFO" "Iniciando Ejercicio 5 — Diseño de Primers"
    log "INFO" "FASTA input: $input_fasta"
    log "INFO" "Config: $config"

    # Validaciones de input
    if [ ! -f "$input_fasta" ]; then
        log "ERROR" "No se encontró el FASTA: $input_fasta"
        exit 1
    fi
    if [ ! -f "$config" ]; then
        log "ERROR" "No se encontró el config: $config"
        exit 1
    fi

    # Validación de formato FASTA básico
    if ! grep -q "^>" "$input_fasta"; then
        log "ERROR" "El archivo $input_fasta no parece ser un FASTA válido (no tiene encabezado >)"
        exit 1
    fi
    log "INFO" "Inputs validados"

    mkdir -p TP_Parte2/Ej5/output

    # Ejecución
    python3 TP_Parte2/Ej5/Ex5.py \
        -i "$input_fasta" \
        -c "$config" \
        -o TP_Parte2/Ej5/output/primers_output.txt >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR" "Ejercicio 5 falló. Revisar log para detalles."
        exit 1
    fi

    # Validación de output
    if [ ! -f "TP_Parte2/Ej5/output/primers_output.txt" ]; then
        log "ERROR" "Output no generado: TP_Parte2/Ej5/output/primers_output.txt"
        exit 1
    fi
    log "OK" "Ejercicio 5 finalizado. Output: TP_Parte2/Ej5/output/primers_output.txt"
}

# =============================================================
# LÓGICA DE SELECCIÓN
# =============================================================

separador
log "INFO" "======= INICIO DE PIPELINE ======="
log "INFO" "Comando: $0 $*"

case "$1" in
    1) ejercicio_1 ;;
    2) ejercicio_2 ;;
    3) ejercicio_3 ;;
    4) ejercicio_4 ;;
    5) ejercicio_5 "$2" "$3" ;;
    all)
        ejercicio_1
        ejercicio_2
        ejercicio_3
        ejercicio_4
        ejercicio_5 "$2" "$3"
        separador
        log "INFO" "======= PIPELINE COMPLETO FINALIZADO ======="
        ;;
    *) usage ;;
esac