![SuppOrth](https://raw.githubusercontent.com/fma91/SuppOrth/main/img/supporth_banner.png)

# SuppOrth

SuppOrth combines orthologue predictions from **OrthoFinder, SonicParanoid2
and Broccoli**. Give it protein FASTAs from two or more species, and it runs
the three programs on the same proteomes, then counts how many programs
identified each cross-species protein pair. Results from a separate **OMA**
analysis can contribute a fourth prediction.

The main output is a table with one row per unique protein pair, the programs
that found it, and alignment statistics where available. One-to-many and
many-to-many relationships are kept as separate pairs.

![SuppOrth workflow](https://raw.githubusercontent.com/fma91/SuppOrth/main/img/supporth_workflow.png)

[Workflow as SVG](https://github.com/fma91/SuppOrth/blob/main/img/supporth_workflow.svg)

## Install

Use Python 3.10–3.12 on macOS or Linux:

```sh
python -m pip install supporth        # Install SuppOrth from PyPI
suppOrth install                     # Install the predictors and their dependencies
suppOrth check                       # Check that the installation is ready
```

The [Python package](https://pypi.org/project/supporth/) provides the SuppOrth
commands. `suppOrth install` then builds the predictors from their source
repositories under `~/.supporth`, with shared copies of their external binaries.
You will need Git and a C/C++ compiler. MAFFT and FastME must already be
available or supplied during installation.

SonicParanoid also needs an indexed Pfam database. Preparing one can take
several hours. If you already have it, pass its location during installation:

```sh
suppOrth install --pfam-db /path/to/profile_db
```

The [installation guide](https://github.com/fma91/SuppOrth/blob/main/supporth_manual.md#installation)
covers Python environments, binary paths and custom installation locations.

## Run

Put one protein FASTA per species or strain in a directory:

```sh
suppOrth run -f proteomes/ -o results/ -t 8
```

Or pass the files directly:

```sh
suppOrth run speciesA.faa speciesB.faa speciesC.faa -o results/ -t 8
```

Protein identifiers must be unique across all input files. SuppOrth uses the
first field of each FASTA header, up to the first space or tab. The predictors
run one after another; `-t 8` requests eight threads for their calculations.

By default, SuppOrth keeps a predicted pair if **at least three tools support
it**, or if its **alignment identity is greater than 80%**. Both thresholds
can be changed:

```sh
suppOrth run -f proteomes/ -o results/ --supportLv 2 --identity_filter 90
```

Support is the number of predictors that found a pair. Use `--supportLv 1`
to keep every predicted cross-species pair.

### Add more strains

To extend a completed analysis, put the original and new FASTAs together and
use the same results directory:

```sh
suppOrth run -f expanded_proteomes/ -o results/ -t 8 --update
```

OrthoFinder and SonicParanoid reuse their earlier sequence searches and
calculate the missing comparisons. Broccoli reruns all four steps. SuppOrth
then rebuilds the consensus for the expanded set of proteomes.

Keep the original FASTAs and predictor settings unchanged, and retain the
`work/` directory containing the previous searches. The
[update guide](https://github.com/fma91/SuppOrth/blob/main/supporth_manual.md#add-proteomes-to-a-completed-run)
explains what is reused and how to include updated OMA results.

## Include OMA

Run OMA separately on the proteomes you want to compare, keeping the same
protein identifiers. Then pass its output to SuppOrth:

```sh
suppOrth run -f proteomes/ -o results/ --oma /path/to/oma_run/Output/
```

SuppOrth reads OMA's `PairwiseOrthologs` tables and includes those pairs in the
support count. OMA installation and execution are handled separately. See
the [OMA guide](https://github.com/fma91/SuppOrth/blob/main/supporth_manual.md#include-oma-predictions)
for accepted folders and help with matching identifiers.

## Results

The main files in `results/` are:

- **`suppOrth.tsv`** contains the retained protein pairs, their support counts,
  a 0/1 column for each tool, and available identity, e-value and bit score.
- **`suppOrth.png` and `suppOrth.svg`** show the overlap between predictors
  before the support and identity filters are applied.
- **`provenance.json`** records the inputs, settings, software versions and
  whether each predictor completed successfully.

Tool labels are **F** (OrthoFinder), **S** (SonicParanoid), **B** (Broccoli) and
**O** (OMA). For example, `SupportedBy=F,S,O` means that OrthoFinder,
SonicParanoid and OMA reported the pair.

![Protein-pair overlaps](https://raw.githubusercontent.com/fma91/SuppOrth/main/img/suppOrth.png)

The [output guide](https://github.com/fma91/SuppOrth/blob/main/supporth_manual.md#output)
includes the first ten rows of an example table and descriptions of the other
output files.

## More options

OrthoFinder uses DIAMOND for similarity searches, MAFFT for multiple sequence
alignment and IQ-TREE for gene trees. Use `--alignment blast` to select BLAST
instead of DIAMOND. You can also pass a predictor's own options in quotes.
For example, to change OrthoFinder's MCL inflation parameter:

```sh
suppOrth run -f proteomes/ -o results/ --opt orthofinder="-I 2.0"
```

Further options, restart commands and installation troubleshooting are in the
[SuppOrth manual](https://github.com/fma91/SuppOrth/blob/main/supporth_manual.md).
