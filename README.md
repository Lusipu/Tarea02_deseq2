# README
## DatasetTitle
Tarea 02: RNA-seq con DESeq2 por consola

### Description
Análisis de expresión diferencial con DESEq2 de la cohorte TCGA-BRCA comparando tejido tumoral ('Primary_Solid_Tumor') frente a tejido normal
adyacente ('Solid_Tissue_Normal') y corrigiendo por la placa de secuencuacion ('plate'). Después se hace enriquecimiento funcional de los genes up y down con GO (Biological Process) y Reactome.

## Organization
**Name:** Universidad de Sevilla

## Researchers
**Name:** Lucía Domínguez García
**Email:** lucdomgar1@alum.us.es
**Subject:** Bioinformática y Análisis Genómico

## General Information

### Keywords
BAG, Bioinformática, Bioquímica

### Language
R

## GitHub Repository
https://github.com/Lusipu/Tarea02_deseq2.git

## Project Structure
```
tarea02_deseq2/
├── README.md
├── environment.yml        # dependencias del entorno conda
├── .gitignore
├── config/
│   └── config.yaml        # rutas, nombres de columnas y umbrales
├── data/
│   ├── raw/               # counts.csv.gz y metadata.csv.gz (aula virtual)
│   └── processed/         # metadata de las muestras usadas (la genera el script)
├── src/
│   └── 01_run_deseq2.R    # script principal
├── results/               # tablas .tsv y figuras .png
└── logs/                  # run_summary.txt y sessionInfo.txt
```

## Instructions

### 1. Datos de entrada

Los datos no están en GitHub (están excluidos en `.gitignore`). Se descargan del aula virtual de la asignatura y hay que colocarlos con estos los siguientes nombres:

```
data/raw/counts.csv.gz
data/raw/metadata.csv.gz
```

`counts.csv.gz` tiene los conteos crudos (genes en filas, muestras en columnas) y `metadata.csv.gz` la información de cada muestra (tipo de tejido y placa).

### 2. Crear el entorno

Hace falta tener conda instalado. En Windows hay que hacerlo dentro de WSL (Ubuntu), porque los paquetes de bioconda no funcionan en Windows directamente.

Desde la carpeta principal del proyecto:

```bash
conda env create -p ./.venv -f environment.yml
conda activate ./.venv
```

El entorno se crea en la carpeta `.venv` dentro del proyecto. La instalación tarda bastante y necesita internet, porque algunos paquetes de anotación (`org.Hs.eg.db`, `GO.db`, `reactome.db`) se descargan desde Bioconductor al instalarse.

Para comprobar que R y DESeq2 funcionan:

```bash
Rscript --version
Rscript -e "library(DESeq2); sessionInfo()"
```

### 3. Ejecutar el análisis

Siempre desde la carpeta principal del proyecto y con el entorno activado:

```bash
Rscript src/01_run_deseq2.R
```

El script va diciendo por qué paso va (del 1 al 8). El paso de DESeq2 es el más lento; con todas las muestras de TCGA-BRCA puede tardar bastante. Si el ordenador tiene varios núcleos se puede poner `nucleos: 4` (por ejemplo) en `config/config.yaml` para que vaya más rápido.

### 4. Resultados

Al terminar se generan en `results/`:

| Fichero | Contenido |
|---|---|
| `deseq2_all_results.tsv` | Todos los genes evaluados, ordenados por padj |
| `deseq2_sig_padj_0.05.tsv` | Genes con padj < 0.05 |
| `deseq2_de_padj_0.05_lfc_1.tsv` | Genes con padj < 0.05 y \|log2FC\| > 1 |
| `top10_up_genes.tsv`, `top10_down_genes.tsv` | Los 10 genes up y down con menor padj |
| `plot_top10_genes.png` | Barras con el log2FC de esos 20 genes |
| `plot_ma.png` | MA plot |
| `plot_volcano.png` | Volcano plot |
| `normalized_counts.tsv` | Conteos normalizados por DESeq2 |
| `top_gene_context.txt` | Nombre completo, log2FC y padj de los genes destacados y de MKI67 |
| `go_bp_up.tsv`, `go_bp_down.tsv` | Enriquecimiento GO BP |
| `reactome_up.tsv`, `reactome_down.tsv` | Enriquecimiento Reactome |
| `plot_go_bp_*.png`, `plot_reactome_*.png` | Dotplots de los 10 términos principales |

Y en `logs/`:

- `run_summary.txt`: número de muestras, genes filtrados, genes DE, términos enriquecidos y duración.
- `sessionInfo.txt`: versión de R y de todos los paquetes.

Para interpretar los resultados: log2FoldChange > 0 significa más expresión en tumor y < 0 más expresión en tejido normal. Se considera diferencialmente expresado un gen con padj < 0.05 (FDR, Benjamini-Hochberg) y |log2FC| > 1, es decir, al menos el doble o la mitad de expresión.

