#!/bin/bash

source ~/bioenv/bin/activate

usage() {
    echo "Uso: $0 {1|2|3|4|5|all} [args...]"
    echo "  1: Procesamiento de Secuencias"
    echo "  2: BLAST"
    echo "  3: Alineamiento Múltiple - MSA"
    echo "  4: Dominios y Motivos PROSITE"
    echo "  5: Diseño de Primers <fasta> [config.json]"
    echo "  all: Corre toda la secuencia (ej5 requiere args extra)"
    exit 1
}

if [ -z "$1" ]; then usage; fi

# --- Definición de pasos ---

ejercicio_1() {
    echo "Ejecutando Ejercicio 1..."
    if [ ! -f "sequence.gb" ]; then
        echo "Error: Falta sequence.gb"
        exit 1
    fi
    python3 ejercicio1.py || { echo "Error en Ejercicio 1"; exit 1; }
}

ejercicio_2() {
    echo "Ejecutando Ejercicio 2 - BLAST..."
    python3 Ejercicio2.py || { echo "Error en Ejercicio 2"; exit 1; }
}

ejercicio_3() {
    echo "Ejecutando Ejercicio 3 - MSA..."
    python3 Ejercicio3.py || { echo "Error en Ejercicio 3"; exit 1; }
}

ejercicio_4() {
    echo "Ejecutando Ejercicio 4 - EMBOSS + PROSITE..."

    if ! command -v patmatmotifs &> /dev/null; then
        echo "Error: EMBOSS no está instalado."
        exit 1
    fi

    if [ ! -f "prosite_data/PROSITE/prosite.lines" ]; then
        echo "Indexando PROSITE por primera vez..."
        mkdir -p prosite_data/PROSITE

        if [ ! -f "prosite_data/PROSITE/prosite.dat" ]; then
            wget https://ftp.expasy.org/databases/prosite/prosite.dat -P prosite_data/PROSITE/
        fi

        if [ ! -f "prosite_data/PROSITE/prosite.doc" ]; then
            wget https://ftp.expasy.org/databases/prosite/prosite.doc -P prosite_data/PROSITE/
        fi

        export EMBOSS_DATA=prosite_data/
        prosextract -prositedir prosite_data/PROSITE/ || { echo "Error en prosextract"; exit 1; }
    else
        export EMBOSS_DATA=prosite_data/
    fi

    python3 Ejercicio4.py || { echo "Error en Ejercicio 4"; exit 1; }
}

ejercicio_5() {
    local input_fasta="${1:-}"
    local config="${2:-primer_config.json}"

    echo "Ejecutando Ejercicio 5 - Diseño de Primers..."

    [[ -z "$input_fasta" ]]   && { echo "Error: Ejercicio 5 requiere un FASTA. Uso: $0 5 <fasta> [config.json]"; exit 1; }
    [[ ! -f "$input_fasta" ]] && { echo "Error: No se encontró $input_fasta"; exit 1; }
    [[ ! -f "$config" ]]      && { echo "Error: No se encontró config: $config"; exit 1; }

    mkdir -p ex5_output

    python3 Ex5.py -i "$input_fasta" -c "$config" -o ex5_output/primers_output.txt \
        || { echo "Error en Ejercicio 5"; exit 1; }
}

# --- Lógica de selección ---

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
        echo "Proceso completo finalizado."
        ;;
    *) usage ;;
esac