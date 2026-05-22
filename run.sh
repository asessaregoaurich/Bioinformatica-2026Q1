#!/bin/bash

# Activar entorno automáticamente
source ~/bioenv/bin/activate

# Función para mostrar que ejercicio se quiere correr
usage() {
    echo "Uso: $0 {1|2|3|4|all}"
    echo "  1: Procesamiento de Secuencias (Ej 1)"
    echo "  2: BLAST (Ej 2)"
    echo "  3: Alineamiento Múltiple - MSA (Ej 3)"
    echo "  4: Dominios y Motivos PROSITE (Ej 4)"
    echo "  all: Corre toda la secuencia"
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

    # Verificar que EMBOSS está instalado
    if ! command -v patmatmotifs &> /dev/null; then
        echo "Error: EMBOSS no está instalado."
        exit 1
    fi

    # Preparar PROSITE si no está indexado
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

# --- Lógica de selección ---

case "$1" in
    1) ejercicio_1 ;;
    2) ejercicio_2 ;;
    3) ejercicio_3 ;;
    4) ejercicio_4 ;;
    all)
        ejercicio_1
        ejercicio_2
        ejercicio_3
        ejercicio_4
        echo "Proceso completo finalizado."
        ;;
    *) usage ;;
esac
