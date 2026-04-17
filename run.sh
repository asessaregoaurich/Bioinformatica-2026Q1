#!/bin/bash

# Activar entorno automáticamente
source ~/bioenv/bin/activate

# Función para mostrar que ejercicio se quiere correr
usage() {
    echo "Uso: $0 {1|2|3|all}"
    echo "  1: Procesamiento de Secuencias (Ej 1)"
    echo "  2: BLAST (Ej 2)"
    echo "  3: Alineamiento Múltiple - MSA (Ej 3)"
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
    python3 ejercicio1.py
}

ejercicio_2() {
    echo "Ejecutando Ejercicio 2..."
    python3 Ejercicio2.py
}

ejercicio_3() {
    echo "Ejecutando Ejercicio 3: MSA..."
    python3 Ejercicio3.py
    echo "Realizando alineamiento múltiple con los mejores hits..."
    python3 Ejercicio3.py || { echo "Error en Ejercicio 3"; exit 1; }
}

# --- Lógica de selección ---

case "$1" in
    1) ejercicio_1 ;;
    2) ejercicio_2 ;;
    3) ejercicio_3 ;;
    all)
        ejercicio_1
        ejercicio_2
        ejercicio_3
        echo "Proceso completo finalizado."
        ;;
    *) usage ;;
esac
