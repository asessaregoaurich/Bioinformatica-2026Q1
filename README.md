# TP Cuatrimestral — Introducción a la Bioinformática

## Enfermedad elegida: Enfermedad de Huntington

Gen: HTT (Huntingtina), cromosoma 4p16.3. Transcripto utilizado: NM_001388492.1 (isoforma 1, MANE Select).

La Enfermedad de Huntington es una enfermedad neurodegenerativa hereditaria autosómica dominante causada por la expansión de repeticiones CAG en el exón 1 del gen HTT. En individuos sanos se observan entre 9 y 35 repeticiones; valores superiores a 40 son patológicos y producen una proteína huntingtina con una región de poliglutamina expandida que es tóxica para las neuronas estriatales.

---

## Estructura del repositorio

```
Bioinformatica-2026Q1/
│
├── inputs/                         ← Archivos de entrada externos (no generados por el pipeline)
│   ├── sequence.gb                 ← GenBank del transcripto NM_001388492.1
│   ├── HTT.fasta                   ← FASTA nucleotídico del transcripto HTT
│   └── primer_config.json          ← Parámetros de diseño de primers (Ej5)
│
├── TP_Parte1/
│   ├── Ej1/
│   │   ├── output/
│   │   │   └── Ej1_ORF.fas         ← 6 marcos de lectura traducidos
│   │   └── ejercicio1.py
│   ├── Ej2/
│   │   ├── output/
│   │   │   ├── blast_local_NM_001388492.1_Frame_+2.xml
│   │   │   └── blast_remote_NM_001388492.1_Frame_+2.xml
│   │   └── Ejercicio2.py
│   └── Ej3/
│       ├── output/
│       │   ├── msa_input.fas
│       │   └── msa_output.aln
│       └── Ejercicio3.py
│
├── TP_Parte2/
│   ├── Ej4/
│   │   ├── output/
│   │   │   ├── Ej4_ORFs.fasta
│   │   │   ├── Ej4_dominios.patmatmotifs
│   │   │   ├── Ej4_pepstats.txt
│   │   │   └── Ej4_resumen.txt
│   │   └── Ejercicio4.py
│   └── Ej5/
│       ├── output/
│       │   └── primers_output.txt
│       └── Ex5.py
│
├── setup.sh                        ← Setup inicial del entorno (correr una sola vez)
├── run.sh                          ← Pipeline principal con logging
├── README.md
└── .gitignore
```

### Dependencias entre ejercicios

| Ejercicio | Input | Fuente |
|-----------|-------|--------|
| Ej1 | `inputs/sequence.gb` | Descargado de NCBI |
| Ej2 | `TP_Parte1/Ej1/output/Ej1_ORF.fas` | Output del Ej1 |
| Ej3 | `TP_Parte1/Ej1/output/Ej1_ORF.fas` + `TP_Parte1/Ej2/output/blast_local_*.xml` | Outputs de Ej1 y Ej2 |
| Ej4 | `inputs/sequence.gb` + `inputs/HTT.fasta` | Descargados de NCBI |
| Ej5 | `inputs/HTT.fasta` + `inputs/primer_config.json` | HTT.fasta generado por Ej4 |

---

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

---

## Ejecución

El script `run.sh` automatiza el flujo completo con validación de inputs/outputs y logging:

```bash
bash run.sh 1        # Solo Ejercicio 1
bash run.sh 2        # Solo Ejercicio 2
bash run.sh 3        # Solo Ejercicio 3
bash run.sh 4        # Solo Ejercicio 4
bash run.sh 5        # Solo Ejercicio 5 (usa inputs por defecto)
bash run.sh all      # Todos en secuencia
```

Cada ejecución genera/actualiza `pipeline.log` en la raíz con timestamps y resultado de cada paso. Este archivo es local y está ignorado por git.

---

## Ejercicio 1 — Procesamiento de secuencias

Script: `TP_Parte1/Ej1/ejercicio1.py`
Input: `inputs/sequence.gb`
Output: `TP_Parte1/Ej1/output/Ej1_ORF.fas`

Lee el mRNA del gen HTT desde el archivo GenBank y genera los 6 marcos de lectura posibles: 3 de la cadena directa (+1, +2, +3) y 3 de la complementaria reversa (-1, -2, -3). Traduce cada uno a aminoácidos y guarda las 6 secuencias en formato FASTA.

El marco de lectura correcto es el **frame +2**. El GenBank indica que la CDS empieza en la posición 146, lo que corresponde a ese frame. Se confirma porque es el único que produce una ORF larga sin codones de stop prematuros, comienza con `MATLEKLMKAFESLKSFQQQQQ...` que coincide con la huntingtina conocida, y el BLAST del Ejercicio 2 lo valida con E-value ~0 e identity 100%.

---

## Ejercicio 2a — BLAST

Script: `TP_Parte1/Ej2/Ejercicio2.py`
Input: `TP_Parte1/Ej1/output/Ej1_ORF.fas`
Output: `TP_Parte1/Ej2/output/`

Puede correr en modo local, remoto o ambos según la configuración. Ejecuta BLASTp contra SwissProt para cada una de las 6 secuencias y guarda un archivo XML por frame en la carpeta `TP_Parte1/Ej2/output/`, con prefijo `blast_local_` o `blast_remote_` según el modo. Los archivos de frames sin hits se eliminan automáticamente — solo se conserva el del Frame +2.

Configuración al inicio del script:

```python
INPUT_FASTA   = "TP_Parte1/Ej1/output/Ej1_ORF.fas"
DB_LOCAL      = "swissprot_db"
OUTPUT_DIR    = "TP_Parte1/Ej2/output"
CORRER_REMOTO = False   # True para correr BLAST remoto contra servidores NCBI
CORRER_LOCAL  = True    # True para correr BLAST local contra SwissProt descargada
```

Para cambiar el modo, editá esas líneas con `nano TP_Parte1/Ej2/Ejercicio2.py` antes de correr. Las combinaciones posibles son:

- `CORRER_LOCAL = True` / `CORRER_REMOTO = False` → corre solo en local, más rápido, requiere SwissProt descargada
- `CORRER_LOCAL = False` / `CORRER_REMOTO = True` → corre solo en remoto contra NCBI, más lento (~2-5 min por secuencia), no requiere DB local
- `CORRER_LOCAL = True` / `CORRER_REMOTO = True` → corre ambos y guarda resultados separados con prefijo `blast_local_` y `blast_remote_`

> Nota: la base de datos `swissprot_db` debe estar en la misma carpeta que `Ejercicio2.py`. No se versiona por su tamaño (~500MB).

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

---

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

---

## Ejercicio 3 — Multiple Sequence Alignment (MSA)

Script: `TP_Parte1/Ej3/Ejercicio3.py`
Input: `TP_Parte1/Ej2/output/blast_local_NM_001388492.1_Frame_+2.xml` + `TP_Parte1/Ej1/output/Ej1_ORF.fas`
Output: `TP_Parte1/Ej3/output/msa_input.fas`, `TP_Parte1/Ej3/output/msa_output.aln`

Requiere MUSCLE instalado.

Descarga las secuencias de los mejores hits del BLAST desde NCBI y las combina con la secuencia query del Frame +2 para construir el archivo de entrada del alineamiento múltiple. El MSA resultante (`msa_output.aln`) incluye la huntingtina humana junto a sus ortólogos en ratón, rata, pez globo y *D. discoideum*.

**Interpretación de Resultados:**
El alineamiento confirma visualmente la alta conservación entre mamíferos y la mayor divergencia con organismos más distantes, consistente con los resultados del BLAST. La secuencia Query (`Human_HTT_Query`) resultó idéntica a la referencia humana (`HD_HUMAN`). Se observa una altísima conservación de la proteína entre los mamíferos (humano, ratón y rata), evidenciada por bloques inmensos de aminoácidos exactamente iguales. Como ejemplo, el extremo N-terminal (que inicia con la secuencia `MATLEKLMKAFESLKSFQQQQ`) se mantiene casi intacto en humanos, roedores e incluso en el pez globo (`HD_TAKRU`).

En contraste, se hace notoria una divergencia extrema en organismos más primitivos como la ameba *Dictyostelium discoideum* (`HD_DICDI`), cuya secuencia presenta numerosos huecos o *gaps* (-) a lo largo de todo el alineamiento. Biológicamente, las extensas regiones que se mantienen idénticas indican la presencia de dominios funcionales críticos que no toleran cambios evolutivos, mientras que las regiones pobladas de *gaps* y variaciones en organismos más distantes demuestran las áreas donde la proteína posee mayor flexibilidad estructural sin perder su función ancestral.

---

## Ejercicio 4 — Dominios y motivos PROSITE

Script: `TP_Parte2/Ej4/Ejercicio4.py`
Input: `inputs/sequence.gb` + `inputs/HTT.fasta`
Output: `TP_Parte2/Ej4/output/`

Lee el archivo GenBank, extrae la secuencia nucleotídica completa y la guarda como FASTA. Luego usa EMBOSS para identificar ORFs (`getorf`), buscar motivos funcionales contra PROSITE (`patmatmotifs`) y calcular estadísticas fisicoquímicas de cada proteína (`pepstats`).

Requiere EMBOSS instalado y PROSITE indexado. El `run.sh` maneja esto automáticamente descargando e indexando PROSITE si no está disponible.

> Nota: los archivos `prosite.dat` y `prosite.doc` (~50MB) se descargan automáticamente y se guardan en `prosite_data/PROSITE/`. No se versionan.

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

---

## Ejercicio 5 — Diseño de primers

Script: `TP_Parte2/Ej5/Ex5.py`
Input: `inputs/HTT.fasta` + `inputs/primer_config.json`
Output: `TP_Parte2/Ej5/output/primers_output.txt`

Diseña primers a partir de la secuencia del transcripto usando una ventana deslizante. Los criterios de diseño están configurados en `primer_config.json`: longitud 18–24 nt, contenido GC 50–60%, Tm máxima 67°C, sin GC en el extremo 3'.

```bash
bash run.sh 5
```

Parámetros de diseño (desde `inputs/primer_config.json`):
- Longitud: 18–24 nt
- Contenido GC: 50–60%
- Tm máxima: 67°C
- Sin GC en extremos 3' ni 5'
- Distancia mínima entre primers: 50 nt
- Cantidad: 5 primers

**Interpretación de Resultados:**

Se obtuvieron 5 primers válidos que cumplen todos los criterios definidos en el archivo de configuración `primer_config.json`. Los cinco resultaron ser de tipo *forward*, con una longitud de 24 nt, un contenido GC del 58.33% y una Tm de 60.8°C. Estos valores homogéneos indican que el espacio de búsqueda en esta secuencia produce candidatos con propiedades termodinámicas muy similares. Aunque la ausencia de primers *reverse* es una limitación para una PCR real (donde se requeriría forzar la selección de al menos un par *forward/reverse*), el algoritmo seleccionado cumplió con priorizar los mejores candidatos individuales basándose estrictamente en las restricciones dadas de Tm y GC.

Una decisión de diseño clave fue excluir la región de repeticiones CAG del exón 1 (posiciones 197–259) mediante el campo `exclude_regions` de la configuración. La consigna requiere que los primers permitan realizar un análisis cuanti y cualitativo de la variante patológica, la cual consiste justamente en la expansión de este tracto de poliglutaminas. Si un primer hibridara dentro de esta región, no sería posible distinguir un alelo sano de uno mutado. Al flanquear la región, el tamaño final del amplicón reflejará el número exacto de repeticiones presentes en el paciente.

Los primers generados respetan este flanqueo estratégico y la separación mínima de 50 nt. Específicamente, los primers #1 (pos. 18–41) y #2 (pos. 171–194) se ubican corriente arriba (*upstream*) del tracto CAG, mientras que los primers #3, #4 y #5 (pos. 674–835) se ubican corriente abajo (*downstream*) de la mutación.

Cabe destacar que los criterios aplicados cubren el diseño básico exigido (longitud, porcentaje GC, Tm máxima y evasión de GC en los extremos). Para un uso experimental in vitro, sería indispensable realizar validaciones adicionales (como la evaluación de formación de *hairpins*, dímeros de primers y un análisis de especificidad mediante BLAST contra el genoma humano para descartar *off-targets*) utilizando herramientas dedicadas como Primer3 o Primer-BLAST.

## Ejercicio 6 - Gen HTT (Huntingtina)

a) Descripción del gen y la proteína
El gen HTT (Gene ID: 3064, cromosoma 4p16.3) codifica la proteína Huntingtina, esencial para el desarrollo neuronal. Una expansión de repeticiones CAG en el exón 1 (>40 repeticiones) genera una región de poliglutamina tóxica que causa la Enfermedad de Huntington, trastorno neurodegenerativo autosómico dominante.

b) Ortólogos
NCBI identificó ortólogos en 810 especies mediante similitud de secuencia. Ensembl identificó 199 especies mediante árboles filogenéticos (más estricto), con 181 relaciones 1 a 1. El gen se conserva exclusivamente en vertebrados (primates, roedores, mamíferos, aves, reptiles y peces), consistente con su función en el sistema nervioso.

c) Transcriptos e isoformas
NCBI reporta 2 transcriptos revisados manualmente (NM_001388492.1 y NM_002111.8). Ensembl identifica 24 transcriptos, de los cuales solo 6 son codificantes de proteína; el resto corresponde a nonsense mediated decay, retained intron o CDS no definida. El transcripto HTT-201 tiene respaldo en ambas bases de datos y corresponde a la isoforma 1 usada en este TP.

d) Interacciones proteína-proteína
La Huntingtina actúa como proteína de andamiaje e interactúa con cientos de proteínas. Sus interactores clave son HAP1/HIP1 (transporte axonal) y caspasas CASP3/CASP6, que al clivar la proteína mutada liberan fragmentos tóxicos. UniProt presenta interacciones curadas experimentalmente; NCBI integra datos de alto rendimiento de BioGRID e IntAct.

e) Gene Ontology 
La Huntingtina se localiza en citoplasma, núcleo y vesículas intracelulares. Participa en transporte axonal, regulación transcripcional, autofagia y desarrollo neuronal. Sus funciones moleculares incluyen unión a proteínas, microtúbulos y factores de transcripción.

f) Vías metabólicas
La HTT es el componente central de la vía de la Enfermedad de Huntington. También interviene en transporte vesicular, autofagia y apoptosis mediada por caspasas.

g) Variantes genéticas 
El SNP rs362307 (3' UTR del gen HTT) está presente en el 65–70% de los pacientes de ascendencia europea, co-segregado con la expansión CAG patogénica. Es utilizado como diana terapéutica para siRNAs que silencian selectivamente el alelo mutado preservando la copia sana.