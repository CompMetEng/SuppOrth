![SuppOrth](https://raw.githubusercontent.com/fma91/SuppOrth/stand_alone/img/supporth_banner.png)

# SuppOrth

SuppOrth runs **OrthoFinder, SonicParanoid2 and Broccoli** on your protein FASTA
files and combines their orthologue predictions. You can also include results
from an existing **OMA** run.

The output contains unique protein pairs between species, the tools supporting
each pair, and alignment statistics when available. It accepts two or more
proteomes and keeps one-to-many and many-to-many relationships.

## Install

Use Python 3.10–3.12 on macOS or Linux. Install SuppOrth, then its predictors:

```sh
python -m pip install supporth        # Install the SuppOrth command-line application
suppOrth install                     # Install the three predictors and prepare their dependencies
suppOrth check                       # Check whether predictors and required binaries are ready
```

SuppOrth is distributed through [PyPI](https://pypi.org/project/supporth/).

The predictors are installed from their source repositories under `~/.supporth`.
You need Git and a compiler toolchain. MAFFT and FastME must already be available
or supplied during installation. SonicParanoid needs an indexed Pfam database;
preparing it can take several hours. To reuse an existing database:

```sh
suppOrth install --pfam-db /path/to/profile_db
```

See the [installation details](https://github.com/fma91/SuppOrth/blob/stand_alone/supporth_manual.md#installation)
for interpreter selection, shared binaries and custom installation locations.

## Run

Put one protein FASTA per species in a directory:

```sh
suppOrth run -f proteomes/ -o results/ -t 8
```

Or pass the files directly:

```sh
suppOrth run speciesA.faa speciesB.faa speciesC.faa -o results/ -t 8
```

Use unique protein identifiers across the FASTAs. SuppOrth reads the first
whitespace-separated field of each header. The predictors run sequentially,
using the requested thread budget.

By default, a pair is kept if **at least three tools support it**, or its
alignment identity is **greater than 80%**. Change the thresholds with:

```sh
suppOrth run -f proteomes/ -o results/ --supportLv 2 --identity_filter 90
```

Use `--supportLv 1` to keep every predicted cross-species pair. Support counts
measure agreement between methods; they are not probabilities.

## Include OMA

Run OMA separately with the same protein identifiers, then supply its output:

```sh
suppOrth run -f proteomes/ -o results/ --oma /path/to/oma_run/Output/
```

SuppOrth imports OMA's `PairwiseOrthologs` tables as a fourth source of support.
It does not install or run OMA. See the [OMA guide](https://github.com/fma91/SuppOrth/blob/stand_alone/supporth_manual.md#include-oma-predictions)
for accepted folders and identifier matching.

## Results

Open these files in `results/`:

- **`suppOrth.tsv`**: protein pairs, support counts, 0/1 tool columns and alignment statistics.
- **`suppOrth.png` and `suppOrth.svg`**: overlap figures from all predictor calls before filtering.
- **`provenance.json`**: inputs, options, software versions and run status.

Tool labels are **F** (OrthoFinder), **S** (SonicParanoid), **B** (Broccoli) and
**O** (OMA). For example, `SupportedBy=F,S,O` means three tools support the pair.

![Protein-pair overlaps](https://raw.githubusercontent.com/fma91/SuppOrth/stand_alone/img/suppOrth.png)

The [SuppOrth manual](https://github.com/fma91/SuppOrth/blob/stand_alone/supporth_manual.md#output)
includes the first ten rows of an example output, all output files, native
predictor options, BLAST/DIAMOND selection and maintenance commands.

## More options

OrthoFinder uses DIAMOND, MAFFT and IQ-TREE by default. Select BLAST with
`--alignment blast`. Pass native predictor arguments in a quoted string:

```sh
suppOrth run -f proteomes/ -o results/ --opt orthofinder="-I 2.0"
```

See the [SuppOrth manual](https://github.com/fma91/SuppOrth/blob/stand_alone/supporth_manual.md)
and [release instructions](https://github.com/fma91/SuppOrth/blob/stand_alone/RELEASING.md).
