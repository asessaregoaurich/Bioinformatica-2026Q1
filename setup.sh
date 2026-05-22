#!/bin/bash

# Detener el script si ocurre un error
set -e

echo "--- Iniciando Setup del Workspace ITBA ---"

#Crear el entorno virtual e instalar dependencias
echo "Configurando entorno Python..."
python3 -m venv ~/bioenv
source ~/bioenv/bin/activate
pip install --upgrade pip
pip install biopython

# Clonar repositorio (si no existe)
REPO_DIR="$HOME/Bioinformatica-2026Q1"
if [ ! -d "$REPO_DIR" ]; then
    echo "Clonando repositorio..."
    git clone https://github.com/asessaregoaurich/Bioinformatica-2026Q1 "$REPO_DIR"
fi
cd "$REPO_DIR"

#Descargar e indexar SwissProt para BLAST local
echo "Preparando base de datos SwissProt..."
if [ ! -f "swissprot" ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
    gunzip swissprot.gz
fi

# Validar si makeblastdb existe antes de ejecutar
if command -v makeblastdb &> /dev/null; then
    makeblastdb -in swissprot -dbtype prot -out swissprot_db
else
    echo "Error: BLAST+ no está instalado. Instálalo para continuar."
    exit 1
fi

echo "--- Setup finalizado con éxito ---"
