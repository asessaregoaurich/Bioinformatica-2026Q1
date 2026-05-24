# TP Cuatrimestral — Introducción a la Bioinformática

## Enfermedad elegida: Enfermedad de Huntington

Gen: HTT (Huntingtina), cromosoma 4p16.3. Transcripto utilizado: NM_001388492.1 (isoforma 1, MANE Select).

La Enfermedad de Huntington es una enfermedad neurodegenerativa hereditaria autosómica dominante causada por la expansión de repeticiones CAG en el exón 1 del gen HTT. En individuos sanos se observan entre 9 y 35 repeticiones; valores superiores a 40 son patológicos y producen una proteína huntingtina con una región de poliglutamina expandida que es tóxica para las neuronas estriatales.

## Setup del entorno

Requiere el Workspace ITBA (BLAST 2.12.0+ ya instalado), Python 3.12+ y conexión a internet. El setup completo puede correrse con un solo script:

```bash
bash setup.sh
```

O paso a paso:

```bash
# Crear entorno virtual e instalar BioPython
python3 -m venv ~/bioenv
source ~/bioenv/bin/activate
pip install biopython

# Clonar el repositorio
cd /home/ITBADC/<tu-usuario>
git clone https://github.com/asessaregoaurich/Bioinformatica-2026Q1
cd Bioinformatica-2026Q1

# Descargar e indexar SwissProt (~137MB comprimido)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
gunzip swissprot.gz
makeblastdb -in swissprot -dbtype prot -out swissprot_db
```

Cada vez que se abre una terminal nueva hay que activar el entorno:

```bash
source ~/bioenv/bin/activate
cd /home/ITBADC/<tu-usuario>/Bioinformatica-2026Q1
```

## Ejecución

El script `run.sh` automatiza el flujo completo:

```bash
bash run.sh 1      # Solo Ejercicio 1
bash run.sh 2      # Solo Ejercicio 2
bash run.sh 3      # Solo Ejercicio 3
bash run.sh all    # Todos en secuencia
```

## Ejercicio 1 — Procesamiento de secuencias

Script: `ejercicio1.py` | Input: `sequence.gb` | Output: `Ej1_ORF.fas`

Lee el mRNA del gen HTT desde el archivo GenBank y genera los 6 marcos de lectura posibles: 3 de la cadena directa (+1, +2, +3) y 3 de la complementaria reversa (-1, -2, -3). Traduce cada uno a aminoácidos y guarda las 6 secuencias en formato FASTA.

```bash
python3 ejercicio1.py
```

El marco de lectura correcto es el **frame +2**. El GenBank indica que la CDS empieza en la posición 146, lo que corresponde a ese frame. Se confirma porque es el único que produce una ORF larga sin codones de stop prematuros, comienza con `MATLEKLMKAFESLKSFQQQQQ...` que coincide con la huntingtina conocida, y el BLAST del Ejercicio 2 lo valida con E-value ~0 e identity 100%.

## Ejercicio 2a — BLAST

Script: `Ejercicio2.py` | Input: `Ej1_ORF.fas` | Output: `blast_results/`

Puede correr en modo local, remoto o ambos según la configuración. Ejecuta BLASTp contra SwissProt para cada una de las 6 secuencias y guarda un archivo XML por frame en la carpeta `blast_results/`, con prefijo `blast_local_` o `blast_remote_` según el modo. Los archivos de frames sin hits se eliminan automáticamente — solo se conserva el del Frame +2.

```bash
python3 Ejercicio2.py
```

Configuración al inicio del script:

```python
INPUT_FASTA   = "Ej1_ORF.fas"
DB_LOCAL      = "swissprot_db"
OUTPUT_DIR    = "blast_results"
CORRER_REMOTO = False   # True para correr BLAST remoto contra servidores NCBI
CORRER_LOCAL  = True    # True para correr BLAST local contra SwissProt descargada
```

Para cambiar el modo, editá esas líneas con `nano Ejercicio2.py` antes de correr. Las combinaciones posibles son:

- `CORRER_LOCAL = True` / `CORRER_REMOTO = False` → corre solo en local, más rápido, requiere SwissProt descargada
- `CORRER_LOCAL = False` / `CORRER_REMOTO = True` → corre solo en remoto contra NCBI, más lento (~2-5 min por secuencia), no requiere DB local
- `CORRER_LOCAL = True` / `CORRER_REMOTO = True` → corre ambos y guarda resultados separados con prefijo `blast_local_` y `blast_remote_`

Resultados obtenidos:

| Frame | Hits | E-value top hit | Identity |
|-------|------|-----------------|----------|
| +1    | 0    | —               | —        |
| +2    | 5    | 0.0             | 100%     |
| +3    | 0    | —               | —        |
| -1    | 0    | —               | —        |
| -2    | 0    | —               | —        |
| -3    | 0    | —               | —        |

Solo el Frame +2 produjo hits significativos, confirmando que es el marco de lectura correcto.

## Ejercicio 2b — Interpretación del resultado BLAST

El BLAST del Frame +2 encontró 5 hits contra SwissProt:

| # | Proteína | Organismo | E-value | Identity |
|---|----------|-----------|---------|----------|
| 1 | Huntingtina (P42858.2) | Homo sapiens | 0.0 | 100% |
| 2 | Huntingtina (P42859.2) | Mus musculus | 0.0 | 91.2% |
| 3 | HD protein homolog (P51111.1) | Rattus norvegicus | 0.0 | 90.8% |
| 4 | HD protein homolog (P51112.1) | Takifugu rubripes | 0.0 | 69.7% |
| 5 | HD protein homolog (Q76P24.1) | Dictyostelium discoideum | 3.19e-18 | 28.8% |

El **score** mide la calidad del alineamiento — un valor de 16594 es altísimo. El **E-value** indica cuántos hits se esperarían por azar en una base de datos de ese tamaño: un E-value de 0.0 significa que la similitud es estadísticamente imposible de explicar por azar, es homología real. El **identity %** indica qué porción de las posiciones alineadas son idénticas.

El hit 1 con 100% de identidad es exactamente la huntingtina humana (UniProt P42858). Los hits 2, 3 y 4 son ortólogos en ratón, rata y pez globo — misma proteína, distinto organismo — lo que refleja alta conservación evolutiva. El hit 5 en *D. discoideum* (un organismo unicelular) con 28.8% de identidad sugiere que esta familia de proteínas tiene un origen evolutivo muy antiguo.

## Ejercicio 3 — Multiple Sequence Alignment (MSA)

Script: `Ejercicio3.py` | Input: `blast_results/blast_local_NM_001388492.1_Frame_+2.xml` + `Ej1_ORF.fas` | Output: `msa_input.fas`, `msa_output.aln`

Descarga las secuencias de los mejores hits del BLAST desde NCBI y las combina con la secuencia query del Frame +2 para construir el archivo de entrada del alineamiento múltiple. El MSA resultante (`msa_output.aln`) incluye la huntingtina humana junto a sus ortólogos en ratón, rata, pez globo y *D. discoideum*.

```bash
python3 Ejercicio3.py
```

El alineamiento confirma visualmente la alta conservación entre mamíferos y la mayor divergencia con organismos más distantes, consistente con los resultados del BLAST.

**[completar interpretación de resultados]**

## Ejercicio 4 — Dominios y motivos PROSITE

Script: `Ejercicio4.py` | Input: `sequence.gb` | Output: `Ej4_ORFs.fasta`, `Ej4_dominios.patmatmotifs`, `Ej4_pepstats.txt`, `Ej4_resumen.txt`

Lee el archivo GenBank, extrae la secuencia nucleotídica completa y la guarda como FASTA. Luego usa EMBOSS para identificar ORFs (`getorf`), buscar motivos funcionales contra PROSITE (`patmatmotifs`) y calcular estadísticas fisicoquímicas de cada proteína (`pepstats`).

```bash
python3 Ejercicio4.py
```

Requiere EMBOSS instalado y PROSITE indexado. El `run.sh` maneja esto automáticamente descargando e indexando PROSITE si no está disponible.

`getorf` identificó 8 ORFs posibles con `-find 1 -minsize 300`. El ORF principal de 3142 aa corresponde a la huntingtina; los 7 restantes son marcos alternativos de menor longitud sin relevancia funcional esperada.

`patmatmotifs` encontró los siguientes motivos en el ORF principal:

| Motivo | Posición | Descripción |
|--------|----------|-------------|
| AMIDATION | 1530–1533 | Posible sitio de amidación C-terminal |
| AMIDATION | 2543–2546 | Posible sitio de amidación C-terminal |
| LEUCINE_ZIPPER | 1444–1465 | Patrón compatible con cremallera de leucina |
| TYR_PHOSPHO_SITE_2 | 2714–2721 | Posible sitio de fosforilación por tirosina quinasa |

Los ORFs 6 y 7 presentaron una coincidencia de amidación cada uno; el resto no mostró hits.

Los sitios de amidación deben interpretarse con cautela: la amidación C-terminal es típica de hormonas peptídicas y proteínas secretadas, no de proteínas intracelulares de gran tamaño como la huntingtina, por lo que probablemente representan coincidencias con el patrón consenso sin relevancia funcional. El motivo leucine zipper también posee baja especificidad y no constituye evidencia concluyente de una región funcional. En cambio, el sitio de fosforilación en tirosina es particularmente interesante dado que la huntingtina sufre múltiples modificaciones postraduccionales que modulan su función, estabilidad y localización celular. En conjunto, estos hallazgos representan evidencia preliminar de posibles regiones funcionales y deben confirmarse con análisis adicionales.

## Ejercicio 5 — Diseño de primers

Script: `Ex5.py` | Config: `primer_config.json` | Input: `HTT_NM_001388492.fasta`

Diseña primers a partir de la secuencia del transcripto usando una ventana deslizante. Los criterios de diseño están configurados en `primer_config.json`: longitud 18–24 nt, contenido GC 50–60%, Tm máxima 67°C, sin GC en el extremo 3'.

```bash
python3 Ex5.py -i HTT_NM_001388492.fasta -c primer_config.json -o primers_output.txt
```

**[completar interpretación de resultados]**
