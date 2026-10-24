![SuppOrth](https://raw.githubusercontent.com/fma91/SuppOrth/main/img/supporth_banner.png)

# SuppOrth manual

SuppOrth runs OrthoFinder, SonicParanoid2 and Broccoli on two or more protein
FASTA files, then combines their predictions into unique cross-species protein
pairs. Results from a separate OMA run can be included as a fourth prediction
set.

Each pair has a support count, one 0/1 column per tool, and alignment statistics
when a matching BLAST or DIAMOND hit is available. If two proteins in one species
are predicted to match one protein in another species, both pairs are retained.
Within-species pairs are excluded.

The predictors provide complementary methods: OrthoFinder uses gene-tree
reconciliation, Broccoli combines graph and phylogenetic inference, SonicParanoid2
uses similarity-based inference with domain analysis, and OMA uses evolutionary
distance consistency. Support counts record agreement between these methods;
they are not calibrated probabilities.

## Contents

- [Installation](#installation)
- [Prepare the FASTAs](#prepare-the-fastas)
- [Run](#run)
- [Include OMA predictions](#include-oma-predictions)
- [Output](#output)
- [Resume and maintenance](#resume-and-maintenance)
- [Troubleshooting](#troubleshooting)
- [Predictor repositories](#predictor-repositories)

## Installation

SuppOrth supports macOS and Linux. Use an existing Python 3.10–3.12
environment and install the application from PyPI:

```sh
python -m pip install supporth  # Install the SuppOrth application
suppOrth install               # Install predictors and prepare their dependencies
suppOrth check                 # Report whether the installation is ready
```

To install the development branch from source instead:

```sh
python -m pip install 'git+https://github.com/fma91/SuppOrth.git@stand_alone'
```

From a local checkout, use `python -m pip install .`. All installation methods
provide the `suppOrth` and `supporth` commands. The Python package does not
bundle predictor repositories, external binaries or the Pfam database;
use `suppOrth install` to prepare them.

`suppOrth install` installs all three predictors sequentially from their upstream
repositories. OrthoFinder is pinned to `v2.5.5`, SonicParanoid2 to `2.0.9`, and
Broccoli to `v1.4.0`. SonicParanoid is compiled from its cloned source using pip.
A compiler toolchain is required for its extensions.

Installations are stored under `~/.supporth`. Set `SUPPORTH_HOME` to use another
location. Helper binaries already on `PATH` are linked into `~/.supporth/bin`;
available upstream releases are downloaded when needed. MAFFT and FastME must
already be installed or supplied with `--adopt`.

```sh
suppOrth install --adopt diamond=/path/to/diamond --adopt mafft=/path/to/mafft
suppOrth install --python /path/to/python3.10
suppOrth install --tools orthofinder --ref orthofinder=v2.5.5
```

OrthoFinder and Broccoli have private virtual environments. SonicParanoid uses a
private package directory with an existing interpreter. A conda-style interpreter
allows this SonicParanoid version to use the shared helper binaries; a different
interpreter may require its exact bundled binary versions. SuppOrth installs
packages with pip and does not install Python or modify the borrowed environment.

SonicParanoid also needs an indexed Pfam database. Installation prepares the
bundled profiles unless an existing database is supplied:

```sh
suppOrth install --pfam-db /path/to/profile_db
```

The database is checked before installation. Its files are hardlinked when
possible, or copied across filesystems. Preparing a new index can take several
hours and requires substantial disk space.

## Prepare the FASTAs

Use one protein FASTA per species. For example:

```text
proteomes/
├── E_coli.faa
├── P_aeruginosa.faa
└── P_putida.faa
```

SuppOrth passes the complete set to each predictor. Adding a third species
therefore changes the inference context; it is not simply a request to repeat
the original two-species comparison. The final table still contains protein
pairs, with species names identifying each comparison.

A header such as `>Ecoli_protein001 description` gives the identifier
`Ecoli_protein001`. Descriptions after the first space are ignored. A header
such as `>sp|P04079|GUAA_ECOLI description` keeps the full identifier
`sp|P04079|GUAA_ECOLI`, including the vertical bars.

If different files contain identifiers such as `protein001`, add species
prefixes before starting. Keep the same identifiers when preparing OMA input,
otherwise its predictions will not match the other results. SuppOrth rejects
duplicate identifiers rather than guessing which species a protein belongs to.

## Run

Pass FASTA files individually or provide a directory:

```sh
suppOrth run speciesA.faa speciesB.faa speciesC.faa -o results/ -t 8
suppOrth run -f proteomes/ -o results/ -t 8
```

Directory input accepts `.faa`, `.fa`, `.fasta`, `.pep` and `.fas`, including gzip
files. Only files directly inside the directory are used. Each predictor runs
once on the complete set of proteomes. The tools run sequentially and manage
their own parallel tasks within the requested thread budget. With no `-t`,
SuppOrth uses the available CPU count minus one.

Protein identifiers are the first whitespace-separated fields in FASTA headers.
They must be unique within each file and across all inputs. Prefix identifiers
with a species name if necessary, and use those same identifiers in any OMA run.
Species labels come from filenames; use `--labels A,B,C` to override them.

For a first run, choose a new output directory and a thread count that leaves
room for other work on your machine. The thread setting controls parallel work
inside each predictor; it does not launch the three predictors together. Larger
proteome sets increase both the number of sequence comparisons and the amount
of downstream grouping and tree inference.

To run only selected predictors, supply their names:

```sh
suppOrth run -f proteomes/ -o results_subset/ --tools orthofinder,broccoli -t 8
```

The support count includes only the tools represented in that run. With two
predictors and no OMA input, no pair can reach a support count of three, so
choose the threshold deliberately.

### Filtering

By default, a pair is retained when at least three tools support it, or when its
alignment identity is strictly greater than 80%.

```sh
suppOrth run -f proteomes/ -o results/ --supportLv 2 --identity_filter 90
```

`--supportLv` accepts 1–4; `--identity_filter` accepts 0–100. A pair without an
alignment hit can still pass the support threshold. To keep every predicted
cross-species pair, use `--supportLv 1`.

The two filter conditions are joined by **OR**. With the defaults, a pair
supported by all three predictors passes even at 50% identity. A pair supported
by one predictor also passes if its identity is 85%. A pair with exactly 80%
identity does not pass the identity condition.

To require support without an identity exception, set the identity threshold
to 100. Since the comparison is strictly greater than the threshold, no valid
percentage identity can satisfy that exception:

```sh
suppOrth run -f proteomes/ -o results_support_only/ --supportLv 3 --identity_filter 100
```

The filter selects among predicted orthologue pairs. A high-identity alignment
alone does not create a new pair that none of the predictors reported.

### Search methods and predictor options

OrthoFinder uses DIAMOND by default, with MAFFT and IQ-TREE for gene-tree
inference. Use `--alignment blast` to select BLAST instead. Its search tables
supply identity, e-value and bit score for the consensus, avoiding an additional
search for those columns.

Broccoli runs its own DIAMOND search because it requires a different output
format. SonicParanoid2 uses its own search workflow. Sharing helper binaries
does not make these predictors share OrthoFinder's search results.

SuppOrth passes the shared DIAMOND and FastTree executable paths directly to
Broccoli. This avoids differences between `FastTree` and `fasttree` filenames
on Linux. Native `-path_diamond` and `-path_fasttree` options take precedence
if you explicitly supply them.

If OrthoFinder is omitted, SuppOrth runs a separate alignment search for the
consensus statistics. `--alignment none` omits these statistics; it does not
switch off the predictors' own searches.

Pass native predictor options as a quoted string. Repeat `--opt` for different
tools:

```sh
suppOrth run -f proteomes/ -o results/ --alignment blast
suppOrth run -f proteomes/ -o results/ --opt orthofinder="-I 2.0"
```

In this example, `-I 2.0` changes OrthoFinder's MCL inflation value. The
surrounding quotes keep the native options together as one SuppOrth argument.
Use the installed predictor's help to check its accepted options and values;
SuppOrth does not translate one tool's parameters into another tool's settings.

Options are passed as arguments, not shell commands. An explicit OrthoFinder
`-S` in `--opt` overrides the search selected by `--alignment`. Method choices
and predictor commits, together with the supplied native arguments, are recorded
in `provenance.json`.

SuppOrth registers threaded search methods in its private OrthoFinder
`config.json`. The search budget is divided among the species comparisons,
including self-comparisons. With `-1`, OrthoFinder searches only one direction.
The original configuration is backed up. Existing
`~/config_orthofinder_user.json` settings remain under the user's control.

## Include OMA predictions

Run OMA separately using the same protein identifiers as the FASTAs supplied to
SuppOrth. OMA is not installed or launched by `suppOrth install` or `suppOrth run`.

Pass the existing OMA run with `--oma`:

```sh
suppOrth run -f proteomes/ -o results/ --oma /path/to/oma_run/Output/
```

The run root or the `PairwiseOrthologs` directory can also be used:

```sh
suppOrth run -f proteomes/ -o results/ --oma /path/to/oma_run/Output/PairwiseOrthologs/
```

A typical OMA output has this layout:

```text
oma_run/
└── Output/
    └── PairwiseOrthologs/
        ├── speciesA-speciesB.txt
        ├── speciesA-speciesC.txt
        └── speciesB-speciesC.txt
```

SuppOrth reads the pairwise tables and matches their protein identifiers against
the input FASTAs. OMA filename abbreviations do not need to match the FASTA
filenames. Predictions involving proteins absent from those FASTAs are excluded.
If OMA identifiers were renamed, restore the corresponding FASTA identifiers
before importing the tables.

An OMA prediction contributes one vote, represented by `O` in `SupportedBy` and
`1` in the `oma` column. A pair reported by OrthoFinder, SonicParanoid and OMA
has `Lv_support=3` and `SupportedBy=F,S,O`. Repeated predictions within one tool
do not add votes. The same filtering rules apply whether three or four tools
are included.

Check the OMA entry in `provenance.json` for the number of imported pairs and
the import status. If no pairs match, check that the protein identifiers agree.
An unreadable OMA result is reported as a failed stage;
other completed predictor results are still written.

## Output

The result directory contains:

- `suppOrth.tsv`: filtered protein pairs, support counts, predictor flags and alignment statistics.
- `suppOrth.png` and `suppOrth.svg`: a Venn diagram, overlap heatmap and exclusive-pair bars from the **unfiltered** predictor calls.
- `agreement_summary.tsv`: counts for each support combination in the filtered table.
- `pairs/`: canonical predictions and JSON dictionaries for each tool.
- `provenance.json`: input checksums, settings, predictor commits, stage status and binary versions.
- `work/`: native predictor outputs and logs.

`SupportedBy` uses comma-separated labels: `F` for OrthoFinder, `S` for
SonicParanoid, `B` for Broccoli and `O` for OMA. The individual tool columns
contain 0 or 1. Alignment fields are empty when no hit was found.

`Lv_support` counts tools, not repeated entries. If one tool reports a pair
several times, it still contributes only one vote. The 0/1 columns use the tool
names (`orthofinder`, `sonicparanoid`, `broccoli`, and `oma` when included), while
`SupportedBy` uses their short labels. Only participating tools have columns.

For example, suppose proteins A1 and A2 both match B1. The output has separate
pairs A1–B1 and A2–B1. It does not add A1–A2, because those proteins belong to
the same species. The ordering of the two proteins in a pair has no biological
direction; A1–B1 and B1–A1 are the same relationship.

`Identity` is percentage identity from a matching search hit. `e_val` and
`Bitscore` describe that alignment. They are not support scores. Empty alignment
fields mean no corresponding hit was available, not zero identity.

The figure and filtered table answer different questions. The figure shows
agreement among the original predictions; the table shows which pairs passed
your thresholds. Their totals can therefore differ.

Example: the header and first ten rows from a two-proteome test:

```text
Species1	Species2	Protein1	Protein2	Lv_support	SupportedBy	orthofinder	sonicparanoid	broccoli	Identity	e_val	Bitscore
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3845011|gb|AAC71645.1|	gi|144227462|gb|AAZ44182.2|	3	F,S,B	1	1	1	53.7	0.0	974.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3844834|gb|AAC71452.1|	gi|144227518|gb|AAZ44306.2|	3	F,S,B	1	1	1	66.9	0.0	960.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|84626139|gb|AAC71565.2|	gi|71852093|gb|AAZ44701.1|	3	F,S,B	1	1	1	38.6	1.17e-301	883.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3844676|gb|AAC71307.1|	gi|144227452|gb|AAZ44165.2|	3	F,S,B	1	1	1	62.8	2.89e-311	858.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3844641|gb|AAC71247.1|	gi|71852011|gb|AAZ44619.1|	3	F,S,B	1	1	1	38.3	7.57e-277	826.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3844923|gb|AAC71566.1|	gi|71852094|gb|AAZ44702.1|	3	F,S,B	1	1	1	37.0	8.2e-271	800.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3844862|gb|AAC71488.1|	gi|71852122|gb|AAZ44730.1|	3	F,S,B	1	1	1	46.2	7.85e-255	722.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|3844728|gb|AAC71355.1|	gi|71851555|gb|AAZ44163.1|	3	F,S,B	1	1	1	55.7	7.06e-238	664.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|1046004|gb|AAC71527.1|	gi|71851549|gb|AAZ44157.1|	3	F,S,B	1	1	1	59.4	4.79e-226	634.0
Mycoplasma_genitalium	Mycoplasma_hyopneumoniae	gi|1045671|gb|AAC71219.1|	gi|71851590|gb|AAZ44198.1|	3	F,S,B	1	1	1	49.7	8.16e-221	625.0
```

![Protein-pair overlaps](https://raw.githubusercontent.com/fma91/SuppOrth/main/img/suppOrth.png)

For two proteomes, `--collapse map1.json,map2.json` also writes
`suppOrth.genes.tsv`. Mapping files can be JSON dictionaries, pickle files or
two-column TSV files. Gene-level results are not plotted.

## Resume and maintenance

### Add proteomes to a completed run

Put the original FASTAs and the new strains in one directory, then use the
existing SuppOrth output directory:

```sh
suppOrth run -f expanded_proteomes/ -o results/ -t 8 --update
```

Include every original proteome with its original label and unchanged file
contents. If you used `--labels`, provide labels for the complete expanded
collection. Repeat the original `--tools`, `--alignment` and `--opt` settings.
Thread counts and consensus thresholds can change. SuppOrth checks input file
hashes, predictor commits and recorded shared-binary versions before starting.

The predictors handle the expanded collection differently:

- **OrthoFinder:** receives its previous working directory with `-b` and only
  the additional FASTAs with `-f`. It reuses previous similarity searches,
  computes the new comparisons and recalculates orthogroups, alignments,
  trees and orthologue predictions. Updated native results are written to a
  separate directory. Earlier working directories must remain available:
  OrthoFinder links to their search files rather than copying them.
- **SonicParanoid:** uses the same native output directory. Its `snapshot.tsv`
  identifies an update; existing alignments and applicable pairwise tables
  are reused. Missing comparisons are calculated and multispecies groups are
  rebuilt. SuppOrth imports only the latest run's final pairwise export.
  Do not pass `-ow`, which forces recomputation. A custom `-p` project ID
  must be new for each update.
- **Broccoli:** reruns all four steps. Its `-steps` option can restart analysis
  from existing intermediate files, but cannot incorporate new proteomes
  into those files. Selective reuse of Broccoli searches is not implemented.

Keep `work/` and the native predictor files; saved `pairs/*.json` files alone
cannot supply the previous searches. If OrthoFinder was omitted, SuppOrth
also reuses its own completed alignment tables between unchanged proteomes.

The update replaces `suppOrth.tsv`, the pair JSONs, summaries and figures with
predictions for the expanded collection. Orthology assignments among the
original species can change. Previous provenance is saved under `history/`;
the new provenance records the added proteomes and native result directories.
Identical search results do not guarantee identical final predictions to a
fresh run: downstream tree inference can produce different trees.
If a predictor fails, the previous consensus files are retained. Inspect its
stage log and retry with `--update`.

If the original run included OMA, supply `--oma` again with an OMA analysis
covering the expanded collection. SuppOrth imports it without running OMA.

### Resume the same inputs

Resume completed stages with the same inputs and settings:

```sh
suppOrth run -f proteomes/ -o results/ --resume
```

Resume reads saved predictions. When provenance is available, changed inputs
are rejected; use `--update` to add proteomes. Resume does not validate changed
predictor options, so keep those unchanged or use a new output directory. Failed stages produce
an error status and a nonzero exit code; completed results remain available.

Inspect helper binaries and installation disk usage:

```sh
suppOrth audit
suppOrth audit --link
suppOrth audit --restore
suppOrth disk
suppOrth disk --clean
```

The audit reports linked, duplicate, divergent and unrunnable binaries.
`--link` keeps original binaries as `.supporth-orig`; `--restore` restores them.
`disk --clean` removes reported build and temporary files. Add `--all` to include
unused predictor checkouts.

## Troubleshooting

### Installation is incomplete

Start with `suppOrth check` and read the failed installation step. Installing
the SuppOrth Python package only supplies the command-line application; the
predictors are prepared by `suppOrth install`. Git, a suitable existing Python
interpreter and a compiler must be available for the source installations.

If a required binary is already installed elsewhere, provide its executable
path with `--adopt`. For a predictor that needs another interpreter, use
`--python TOOL=/path/to/python`. SuppOrth does not install Python itself.

### SonicParanoid is preparing Pfam for a long time

Building the profile index is a substantial operation. If an indexed database
already exists, supply its `profile_db` directory with `--pfam-db`. Copying only
one index file is not enough: the associated MMseqs database files must also
be present. Allow space for a copy if the source and destination are on
separate filesystems.

### A run finishes with an error

Check the stage status in `provenance.json` and the logs under `work/`. SuppOrth
keeps results from completed stages, so an output table can exist even when
another predictor failed. The nonzero exit status identifies an incomplete run;
inspect the failed stage before interpreting its support counts.

A tree-building step can be slow without being stuck. Check its native logs
and process activity before interrupting it. Runtime depends on sequence
content as well as the number of proteins.

### OMA contributes no pairs

Compare the protein identifiers in its pairwise tables with the first fields
of your FASTA headers. Similar-looking accessions are not necessarily identical
identifiers: `P04079` and `sp|P04079|GUAA_ECOLI` are different strings. Also check
that the supplied folder contains the pairwise output for the same proteomes.

### The filtered table has fewer pairs than expected

Check the support threshold, identity threshold and available alignment hits.
A support level of four requires OMA as well as all three installed predictors.
To inspect all reported pairs, use `--supportLv 1` in a new output directory.
Keep the figure's unfiltered counts separate from the filtered table's totals.

## Predictor repositories

- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- [SonicParanoid2](https://gitlab.com/salvo981/sonicparanoid2)
- [Broccoli](https://github.com/rderelle/Broccoli)
- [OMA standalone](https://omabrowser.org/standalone/)
