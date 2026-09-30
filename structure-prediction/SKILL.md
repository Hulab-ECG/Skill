---
name: structure-prediction
description: Retrieve, model, or predict protein structures for candidate protein sequences using template search against PDB or local AlphaFold DB sequence search, homology modeling, or user-specified predictors, with pLDDT confidence filtering and resumable batch processing.
---

# Structure Prediction

Given protein sequences, retrieve known structures first and predict missing
structures when required.

The default workflow is resumable and provenance-first. For large bacterial
genomes or pan-proteomes, use template search and AFDB retrieval to cover most
sequences before considering expensive de novo prediction.

## 1. Run the workflow

Core scripts live at `scripts/` in this skill directory.

### Option A: use the wrapper script

Export the standard input paths, then run:

```bash
export STRUCTURE_PRED_ROOT=/path/to/workdir
export STRUCTURE_PRED_AFDB_FASTA=/path/to/afdb_sequences.fasta

bash <skill_directory>/scripts/run_structure_pipeline.sh \
  --pdb-db  /path/to/pdb.dmnd \
  --afdb-db /path/to/afdb.dmnd \
  --threads 44 \
  --model-source pdb
```

Direct AFDB models and AFDB homology templates must pass the AFDB quality gate
(`mean pLDDT >= 60`). The default threshold can be overridden only with:

```bash
export STRUCTURE_PRED_AFDB_MIN_PLDDT=60.0
```

`--model-source` controls homology modeling:

- `pdb` runs only PDB homology modeling after direct downloads (default).
- `afdb` runs only AFDB homology modeling.
- `all` runs both PDB and AFDB homology modeling sequentially.
- `none` stops after direct PDB/AFDB structure materialization.

### Option B: run each stage manually

```bash
export STRUCTURE_PRED_ROOT=/path/to/workdir
export STRUCTURE_PRED_AFDB_FASTA=/path/to/afdb_sequences.fasta

# The HOMD standard layout is expected under $STRUCTURE_PRED_ROOT:
#   queries/all.faa
#   queries/query_map.tsv
#   paper_diamond/pdb_q50_s50/pdb_hits_q50_s50.tsv
#   paper_diamond/afdb_q50_s50_v2/afdb_hits_q50_s50_v2.tsv

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# 1. Build the query manifest (optional if queries/all.faa already exists)
python3 <skill_directory>/scripts/build_query_manifest.py

# 2. PDB template search
diamond blastp --ultra-sensitive --query-cover 50 --subject-cover 50 --id 50 --iterate \
  --query "$STRUCTURE_PRED_ROOT/queries/all.faa" --db /path/to/pdb.dmnd \
  --out "$STRUCTURE_PRED_ROOT/paper_diamond/pdb_q50_s50/pdb_hits_q50_s50.tsv" \
  --outfmt "$OUTFMT" --threads 44 --max-target-seqs 200

# 3. Classify PDB hits into direct and homology candidates
python3 <skill_directory>/scripts/analyze_pdb_hits.py

# 4. Download direct PDB structures first
python3 <skill_directory>/scripts/download_pdb_direct.py

# 5. AFDB template search
diamond blastp --ultra-sensitive --query-cover 50 --subject-cover 50 --id 50 --iterate \
  --query "$STRUCTURE_PRED_ROOT/queries/all.faa" --db /path/to/afdb.dmnd \
  --out "$STRUCTURE_PRED_ROOT/paper_diamond/afdb_q50_s50_v2/afdb_hits_q50_s50_v2.tsv" \
  --outfmt "$OUTFMT" --threads 44 --max-target-seqs 200

# 6. Classify AFDB hits into direct and homology candidates
python3 <skill_directory>/scripts/analyze_afdb_hits.py

# 7. Download direct AFDB structures first
python3 <skill_directory>/scripts/download_afdb_direct.py

# 8. Build the homology-modeling pending manifest (merged, PDB-only, AFDB-only)
python3 <skill_directory>/scripts/build_homology_pending.py

# 9. Run homology modeling after direct downloads complete
python3 <skill_directory>/scripts/run_homology_modeling.py \
  --manifest "$STRUCTURE_PRED_ROOT/homology_new/homology_pending_pdb_only.tsv" --threads 6

# 10. Optional retry for AFDB DOWNLOAD_FAILED lists
python3 <skill_directory>/scripts/retry_afdb_direct.py
```

### Core scripts

| Script | Purpose |
| --- | --- |
| `scripts/build_query_manifest.py` | Build `queries/all.faa` and `queries/query_map.tsv` from `manifests/recognized.faa.list` |
| `scripts/analyze_diamond_hits.py` | Generic classifier for 14-column DIAMOND hits into direct and homology candidates |
| `scripts/analyze_pdb_hits.py` | Classify PDB DIAMOND hits into direct and homology candidates |
| `scripts/analyze_afdb_hits.py` | Classify AFDB DIAMOND hits into direct and homology candidates |
| `scripts/download_pdb_direct.py` | Download direct PDB chains, falling back from PDB format to mmCIF conversion |
| `scripts/download_afdb_direct.py` | Download direct AFDB models from API `pdbUrl`, with gzip and pLDDT support |
| `scripts/afdb_utils.py` | Shared AFDB mean-pLDDT parser and quality gate |
| `scripts/audit_afdb_structure_quality.py` | Audit/quarantine direct AFDB outputs below the pLDDT threshold |
| `scripts/retry_afdb_direct.py` | Retry previously failed AFDB direct downloads |
| `scripts/build_homology_pending.py` | Build merged, PDB-only, and AFDB-only homology manifests excluding direct hits |
| `scripts/filter_homology_manifest.py` | Filter a merged homology manifest by `PDB` or `AFDB` source |
| `scripts/extend_diamond_hits.py` | Append two columns, `qlen`/`slen`, turning older 12-column DIAMOND hits into 14-column hits |
| `scripts/run_homology_modeling.py` | Run MODELLER homology modeling from a candidate manifest |
| `scripts/run_structure_pipeline.sh` | Wrapper that chains DIAMOND search, analysis, direct download, and homology stages |

Paths are read from the `STRUCTURE_PRED_ROOT` and
`STRUCTURE_PRED_AFDB_FASTA` environment variables. The scripts do not depend
on any fixed server user, host, or database path, so they can be reused on
another machine by pointing these variables to your local work directory and
AFDB FASTA snapshot.


## Input and Output Conventions


- Input is one or more FASTA protein-sequence files, normally one file per
  species/genome.
- Default input rule: include every non-empty `.faa` file; exclude only files
  matching `*.hypotheticals.faa` (Bakta-style hypothetical protein files).
- Treat aliases such as `GCA_*.1.faa`, `GCF_*.1.faa`, `GCA_*.1_protein.faa`,
  and `GCF_*.1_protein.faa` as protein FASTA files. Do not infer inclusion from
  the file-name suffix alone; validate that the file contains protein FASTA and
  is non-empty.
- For every protein, write at most one final structure file:
  `structure/<protein_sequence_id>.pdb`, placed in the same species directory
  as its source FASTA file.
- `protein_sequence_id` is the first whitespace-separated field of the FASTA
  header, after the leading `>`.
- If a FASTA file has duplicate sequence IDs, resolve or flag before starting.
- The same `protein_sequence_id` must be used for PDB, AFDB, and Phase B
  outputs so downstream coverage checks are unambiguous.
- Always finish with a per-file/per-species summary that reports how many
  proteins have known structures, predicted structures, and no structure.
- Preserve the original FASTA identifier in all mapping, hit, candidate, and
  provenance tables. Use a stable internal query key only when external IDs are
  not unique, but always retain the original mapping.
- Write temporary files atomically, then move them into the final
  `structure/` path. Never leave an empty or partially written `.pdb` file in
  the final location.
- Every final structure must have provenance: query ID, source FASTA, template
  or model ID, source type, method, identity/coverage when relevant, quality
  metric, download or run timestamp, and output directory.

## Resumability and Resource Limits

- Maintain a global query FASTA and query map before starting large searches.
- Make each phase independent and idempotent. Before running, skip proteins
  with an existing valid, non-empty final structure file.
- Append structured status to a phase log rather than rewriting it during the
  run. Use one status per query or one per API result, and retain all retries.
- Do not mark a phase complete until every query is accounted for. Record
  `FOUND`, `EXISTS`, `MISSING`, `LOW_QUALITY`, `NO_TEMPLATE`, `PENDING`,
  `USER_LIMIT`, or an explicit `ERROR` with the reason.
- For long runs, write a periodic status/counter file. Display processed,
  remaining, succeeded, low-quality, and error counts without locking jobs.
- Let the user set thread counts. Before using more than a modest number of
  threads, confirm the server load limit; never assume exclusive access.
- Use bounded workers and exponential backoff for remote APIs. Rate-limit
  requests rather than restarting from scratch after bursts of 404s or errors.
- Separate permanent nonexistence from transient failures. A 404 is usually a
  missing entry; timeouts, 429, 5xx, and incomplete responses require retry.
- Resume from logs and existing files after interruption; do not silently
  restart completed searches.

## Workflow Overview

Use eight ordered stages grouped as:

- **Phase A — Known/template structures**: Stages 1-6.
- **Phase B — De novo prediction**: Stage 7.
- **Phase C — Final selection and reporting**: Stage 8.

The final structure priority remains: experimental PDB > direct AFDB model
passing `mean pLDDT >= 60` > homology model > de novo prediction.

### Stage 1. Prepare a Protein Manifest

1. Discover and validate all qualifying `.faa` files.
2. Collapse duplicate protein IDs within each source FASTA or mark them as a
   preflight error.
3. Build a normalized query FASTA, query-to-file map, and sequence checksum
   table. All downstream outputs must reference these original IDs.
4. Establish an existing-structure inventory and define `pending` as the set
   without a valid structure file.

### Stage 2. Structure Sequence Databases

Let the user choose one of two database sources:

1. **Existing local library**: accept an explicit `--pdb-db-dir`
   (or equivalent configuration) and locate either a ready-to-use DIAMOND
   database (`*.dmnd`) or a PDB seqres FASTA
   (`pdb_seqres.txt` / `pdb_seqres.txt.gz`).
2. **Download and build locally**: when no local library is supplied and the
   user permits downloads, retrieve the current RCSB PDB seqres FASTA into a
   project-specific database directory, decompress it if needed, and build a
   DIAMOND database with `diamond makedb`.

Official source reminder:

```text
https://files.rcsb.org/pub/pdb/derived_data
```

`pdb_seqres.txt.gz` is normally under this RCSB directory. After download,
decompress it before `diamond makedb`.

Rules:

- Never invent or overwrite a database directory. If the path does not exist,
  ask for the intended path or obtain explicit permission to create it.
- Validate the selected library before searching: required files exist and are
  non-empty, the FASTA contains PDB accession headers, and a DIAMOND database
  is readable by the installed DIAMOND version.
- If a FASTA library exists without a matching DIAMOND database, build the
  DIAMOND file in the same user-specified database directory.
- Record the database directory, input FASTA path, DIAMOND database path,
  checksums, sequence count, PDB IDs present, snapshot/download date, and
  DIAMOND version in the run provenance.
- Rebuild only when the library is missing, unreadable, obviously stale, or the
  user requests a refresh. Use a filesystem lock for long builds so concurrent
  runs cannot corrupt the database.
- Use user-controlled thread and disk limits; PDB seqres expands significantly
  while DIAMOND is being created.

#### AFDB

The primary AFDB route is a local AFDB full-sequence FASTA searched with
DIAMOND, for example `--afdb-db-dir /path/to/afdb`. The
UniProt-mapping / AFDB-API route is a fallback, not the default path.

Official source reminder:

```text
https://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta
```

1. Accept an explicit AFDB database directory containing `sequences.fasta` or
   another AFDB full-sequence FASTA, and preferably a `.dmnd` database.
2. If only FASTA exists, build the DIAMOND database with `diamond makedb` in
   the same directory after the user has approved the resource usage.
3. Record the FASTA and DIAMOND paths, snapshot/version, file sizes, checksums,
   sequence count, and DIAMOND version.
4. Do not index an incomplete download. Validate the FASTA terminator and
   sequence count before building the index.

Snapshot integrity rules:

- Record both file size and download date, but do not use file size alone as a
  completeness check. Validate the FASTA terminator, search for malformed or
  empty records, and verify the reported sequence count before building the
  DIAMOND index.
- Label every search result with the local AFDB snapshot path and snapshot date.
  Do not merge reads that came from a local snapshot with fallback API results
  unless provenance makes the snapshot/API distinction explicit.
- If API lookup fails while the local AFDB FASTA has a matching accession,
  report `LOCAL_SNAPSHOT_PRESENT_API_404` and use the locally indexed evidence
  after structural materialization succeeds.
- Do not guess AFDB file versions from common patterns such as
  `model_v1`, `model_v4`, or `model_v6`. The file URL differs by entry and
  release; always fetch the official API response once to obtain the exact
  `entryId` and `pdbUrl` before downloading.
- Some AFDB entries found on the website or in the local AFDB FASTA use
  numeric identifiers (for example `AF-0000000008184908`) that are accepted by
  the API as a full entry ID but not as a UniProt-style `AF-...-F1` accession.
  Treat these as valid AFDB lookup keys even when they do not follow the
  standard UniProt accession pattern.
- AFDB `pdbUrl` payloads may be gzip-compressed. Check for the gzip magic bytes
  at runtime, decompress before validating the PDB header, and accept either
  `HEADER`, `ATOM`, `CRYST1`, or `MODEL` at the start of a valid PDB document.

### Stage 3. Template Search Against PDB

Stage 3 is not complete after the DIAMOND search or intermediate assignment. It
is only complete after every direct-reuse PDB structure has been materialized as
`structure/<protein_sequence_id>.pdb` in the source species directory.

#### 3a. Run template search

Use a configurable sequence-search strategy:

- Preferred large-scale local search: DIAMOND `blastp` against a local PDB
  seqres database supplied by or built for the user.
- Accept `--pdb-db-dir` as the first database argument for the run.
- If a local PDB database is not available, use a remote sequence search API
  before falling back to identifier-based UniProt cross-reference retrieval.
- Canonical search parameters:
  `diamond blastp --ultra-sensitive --query-cover 50 --subject-cover 50
  --id 50 --iterate`.
- Store all threshold-passing hits, not just the first hit, so selection can be
  rerun without repeating the search.

#### 3b. Classify hits

- Direct reuse: full-length 100% identity and 100% bidirectional coverage.
- Template-only: passes search thresholds but is not exact/full length. Used
  later by homology modeling, not final structure reuse.

#### 3c. Materialize direct-reuse structures (MUST NOT SKIP)

For every direct-reuse PDB chain:

1. Parse the subject ID, e.g. `1abc_A` gives `pdb_id=1abc`, `chain=A`.
2. Download `https://files.rcsb.org/download/{PDB_ID}.pdb` first.
3. If legacy `.pdb` returns 404, retry
   `https://files.rcsb.org/download/{PDB_ID}.cif`.
4. Parse `.pdb` or `.cif` with a mature parser:
   - `Bio.PDB.PDBParser` / `Bio.PDB.MMCIFParser`, `gemmi`, or equivalent.
   - Do not text-munge mmCIF coordinates by hand.
5. Extract the selected chain. For large complexes or assemblies, do NOT copy
   the whole complex unless the matched target is the whole entry.
6. Handle multi-character author chain IDs by mapping to single-character PDB
   chain IDs using the parser.
7. Write the final file atomically:
   - target: `<species_dir>/structure/<original_protein_sequence_id>.pdb`
   - write temp file in the same directory, then `os.replace`.
   - never leave an empty or partial `.pdb` in the final path.
8. Record provenance in `pdb_direct_download_status.tsv`.

Do not advance to Stage 4 until every exact/full-coverage hit has been attempted
and accounted for in the status log.

- DIAMOND is only template discovery. Final assignment is made in the template
  evaluation step. Keep DIAMOND identity and coverage columns separate from
  final selection metrics.


### Stage 4. AFDB Template Search and Candidate Retrieval

The primary route mirrors PDB template search:

1. Search the local AFDB full-sequence DIAMOND database for every protein still
   missing a structure. Use the same configurable search parameters as Stage 3
   unless the user explicitly supplies different thresholds.
2. Save all threshold-passing AFDB hits, not only the highest-scoring hit.
3. Classify a hit as direct AFDB reuse only when it has 100% sequence identity,
   full bidirectional coverage, and mean pLDDT `>= 60`. Materialize a passed
   model under the original `protein_sequence_id`.
4. Reject direct AFDB models that do not meet the pLDDT quality gate.
5. Treat any non-direct threshold-passing hit as template-only and route it to
   homology-modeling candidate evaluation; do not use it as a direct AFDB
   structure.

The fallback route is mapping/API retrieval:

1. Map original proteins to candidate UniProt accessions while preserving all
   mapping sources and scores.
2. For protein IDs not accepted by a RefSeq idmapping channel, recover xrefs
   from genome annotation or bridge through UniParc/sequence checksums.
3. Query AFDB by each candidate UniProt accession. Query all candidates when
   multiple mappings exist, then select one best model for the original query.
   AFDB is accession-indexed, so do not infer sequence equivalence without
   mapping provenance or identity/coverage evidence.
4. AFDB can complement, but not override, an already assigned valid PDB
   structure. Only use it for queries still missing a structure unless the user
   explicitly asks for a re-selection pass.

### Stage 5. PDB/AFDB Intermediate Assignment and Routing

Evaluate PDB and AFDB candidates together, but do not yet make the full final
assignment. This step only assigns a direct highest-priority result when one
exists and routes everything else to the next stage.

For each original protein ID:

1. **Experimental PDB**: assign the best valid passing PDB chain if available.
2. **AFDB direct model**: assign the best direct AFDB model only when no valid
   PDB result was assigned and the model passes `mean pLDDT >= 60`.
3. **Homology candidate**: if neither direct source produced a structure but a
   template passes the homology-modeling identity/coverage thresholds, mark it
   for Stage 6 homology modeling.
4. **Phase B pending**: if no direct structure exists and no qualifying
   template is available, mark it for Stage 7.

Intermediate requirements:

- Resolve one-to-many mappings by considering all candidates and keeping one
  best output per original protein sequence ID.
- Never name a final file by UniProt, PDB, or AFDB accession; convert it back to
  the original `protein_sequence_id`.
- When several candidates pass, choose deterministically. Recommended order:
  exact identity, higher bidirectional sequence coverage, higher UniProt/AFDB
  quality, better PDB method/resolution, fewer template gaps, stable accession,
  then earlier report order.
- Write exact, complete, and unambiguous provenance for the selected candidate.

Do not run de novo prediction in Stage 5. A protein must first pass through
the homology-modeling eligibility rule.



## MODELLER Environment and Configuration

MODELLER is the default homology-modeling engine in this skill when the
user's pipeline accepts it and a license is available.

Example installation (replace every path with your local MODELLER):

- Conda environment: `meta` or the environment with MODELLER installed.
- Version: for example `MODELLER 10.8`.
- Executable: `mod10.8` inside the active environment.
- Python import: `conda run -n meta python -c "import modeller; print(modeller.info.version)"`.
- Modeller config file:
  `$CONDA_PREFIX/lib/modeller-10.8/modlib/modeller/config.py` when running
  through the active Conda environment.
- Configured license key: keep it in the local MODELLER `config.py`, not in
  version control.

Run MODELLER through the configured Conda environment:

```bash
# Or simply activate the environment and run `mod10.8 build_model.py`.
conda run -n meta python build_model.py
```

`MODELLER 10.8` can read mmCIF templates directly. Use `.pdb` when the source
provides it and `.cif` is acceptable when the template only exists as mmCIF.

License/configuration notes:

- Do not commit MODELLER license keys to version control.
- Before launching homology modeling, check that `import modeller` succeeds and
  that the configured license key returned by Modeller is valid.
- If the license value in `config.py` is still `XXXX`, stop and ask for the
  correct key before running models.
- Re-run the import smoke test after any MODELLER upgrade or environment
  rebuild.

### Stage 6. Homology Modeling

1. Build the homology-modeling pending list: `identity > 90%`, query coverage
   `>= 90%`, and subject coverage `>= 90%`.
2. Run homology modeling only if the user has enabled it and the required
   dependencies/license are available.
3. Materialize and validate the selected template chain before invoking the
   modeler. A DIAMOND hit is not usable evidence by itself; the template file
   and chain must be present locally.
4. Keep every successful output and its provenance recorded under its original
   protein ID.
5. Re-evaluate the provisional Stage 7 pending list after homology modeling:
   proteins with a passing homology model are removed from Stage 7, while
   proteins without one remain there according to the user's fallback policy.

Template handling rules:

- Cache retrieved templates in a run-level template directory, for example
  `homology_templates_pdb/{pdb_id}_{chain}.pdb` or
  `homology_templates_afdb/{accession}.cif`. Write cache entries atomically and
  validate them before reuse. Never trust a zero-byte or truncated cache file.
- For a PDB template, try the legacy PDB coordinate file first. If it returns
  HTTP 404, immediately retry the mmCIF document. Many newer PDB entries are
  mmCIF-only; a `.pdb` 404 is therefore not by itself a permanent template
  failure.
- MODELLER 10.8 and newer can read mmCIF directly. Do not reject a template only
  because the original file is `.cif`. For older modelers, convert mmCIF with a
  mature tool such as `gemmi` before modeling.
- Parse the template with a mature structural parser and verify that the
  selected chain and expected residues are present. Record whether the working
  template was `.pdb` or `.cif`.
- Classify template failures explicitly:
  `TEMPLATE_404` after both coordinate formats have been tried,
  `TEMPLATE_TRANSIENT` for DNS/timeout/HTTP 429/HTTP 5xx,
  `TEMPLATE_INVALID` for malformed, empty, wrong-chain, or sequence-incompatible
  templates, and `MODEL_FAILED` for modeler execution errors. Retry only
  transient failures unless the user asks for a different retry policy.
- Route a failed template/model to Phase B only after retry policy has finished
  and the user's fallback policy permits it.

### Stage 7. Phase B: De Novo Structure Prediction

1. Build the final Stage 7 (Phase B) pending list after homology modeling:
   proteins with no valid PDB, AFDB, or homology-model output.
2. Run the user-specified Phase B predictor for these proteins.
3. Keep every successful output and its provenance recorded under its original
   protein ID.
4. Record low-confidence, failed, and user-limited runs separately so they can
   be resumed without reevaluating completed proteins.

### Stage 8. Unified Candidate Selection and Reporting

After PDB retrieval, AFDB retrieval, homology modeling, and de novo prediction
are complete, first select at most one final structure per original protein ID
in this exact priority order:

1. **Experimental PDB**
2. **AFDB direct model** passing the AFDB mean-pLDDT gate
3. **Homology model**
4. **Phase B predicted structure**

Rules:

- A valid output assigned in Stage 5 is normally retained; rerun final
  selection only after a retrieval or structural failure, a user-requested
  refresh, or a provenance inconsistency.
- If multiple candidates pass, use the deterministic ordering above, then
  exact identity, higher bidirectional sequence coverage, higher source quality,
  better PDB method/resolution, fewer template gaps, stable accession, and
  earlier report order.
- Save the selected file as `structure/<protein_sequence_id>.pdb` and record
  why it was chosen.
- Cancel or rename superseded trial outputs only when the final selection is
  complete; do not leave empty or partial structure files.

After final selection, generate coverage and provenance reports. Reporting must
not change structure assignment.

1. Build a per-FASTA, per-species, and global coverage table.
2. Report counts by final source: PDB, AFDB, homology model, Phase B, missing,
   and pending.
3. Include low-confidence, failed, no-template, dependency-unavailable, and
   user-limited cases separately.
4. Output a machine-readable TSV in addition to the human-readable summary.
5. Regenerate the report from provenance logs and existing final structure
   files without rerunning retrieval.

## Phase A Implementation Notes

### Known-Structure Retrieval

Purpose:

- Retrieve known structures before running new predictions.

Priority:

1. PDB experimental structures.
2. AlphaFold DB predicted structures with `mean pLDDT >= 60`.
3. Qualified homology-modeling templates.
4. Phase B prediction.

Feeding rule: do **not** send proteins that already obtained a PDB structure to
AlphaFold DB. Do not send proteins with a known PDB, AFDB, homology-model, or
Phase B structure to another prediction stage.

Steps:

1. Build the query manifest and wait for the preflight check to pass.
2. Run template search against PDB.
3. Run AFDB retrieval for proteins that still need structures.
4. Perform PDB/AFDB intermediate assignment and routing.
5. Run homology modeling for proteins that qualify.
6. Run the user-specified Phase B predictor for remaining proteins.
7. Perform unified candidate selection and final assignment.
8. Generate per-species and global coverage reports.

If the user asks for the paper-style workflow, use these exact thresholds:

- Exact direct reuse: 100% sequence identity and full bidirectional coverage.
- Homology modeling: identity greater than 90%, query coverage at least 90%,
  and subject coverage at least 90%.
- Do not model or predict proteins that already received a valid structure.

For papers, the homology-modeling criterion is often written as
`identity > 90%` with `90% coverage`, but "coverage" must be made explicit as
query coverage, subject coverage, or bidirectional coverage.

### UniProt mapping notes

- `RefSeq_Protein` accession mapping is straightforward for WP_/NP_/YP_ style
  IDs.
- INSDC-style IDs such as `CAP`, `MFC`, and `CAQ` are often **not** accepted by
  the `RefSeq_Protein` idmapping channel. In that case, recover UniProt
  candidates from the genome annotation (`gbff`/`gff3`) xrefs, e.g. UniParc
  `UPI...` or UniRef entries, or bridge through sequence checksums/UniParc
  searches.
- Bakta annotations may include `RefSeq:WP_...`, `UniParc:UPI...`,
  `UniRef:UniRef100_...`, `UniRef:UniRef50_...`, and
  `UniRef:UniRef90_...` xrefs. Use RefSeq where available for direct ID
  mapping, and UniParc as a bridge when direct mapping fails.
- UniParc hits can return multiple UniProt entries. Store all candidates and
  let AFDB evaluation choose the best passing model; do not discard ambiguity
  during mapping.
- Merge UniParc-derived candidates back into the canonical UniProt mapping
  table rather than treating them as a separate workflow output.
- UniProt idmapping stream jobs do not support extra `fields` parameters in the
  batch submit request; adding `?fields=` can return an empty `{results: []}`.
- The old `/idmapping/status/{jobId}` endpoint may return 404; use the current
  job status/result polling workflow and allow long wait times for large batches.
- For very large inputs, split into chunks, write mapped chunks separately, and
  merge later.
- Some UniProt candidates return AFDB 404 because the accession has no AFDB
  model. This is a valid negative result, not a fatal search error. It should
  be distinguishable from transport, rate-limit, malformed-response, and
  download failures. If later API releases add or remove entries, rerun the
  mapping/API pass only under a clearly labeled new AFDB release; do not alter
  provenance from the old run.

### PDB retrieval notes

- Preferred PDB method ordering: X-ray/crystallography > Cryo-EM > NMR > other.
- Within the same method, prefer lower resolution.
- If a UniProt entry maps to multiple PDB IDs, keep only one best PDB model per
  protein and save it as `structure/<protein_sequence_id>.pdb`.
- Record method, resolution, chains, and source PDB ID in a retrieval log.
- Keep mmCIF as a stable cache or working document when useful, but convert the
  selected compatible chain to PDB when producing the final structure because
  the contract is `structure/<protein_sequence_id>.pdb`.
- When converting mmCIF to PDB, handle long structure sizes atomically and
  clamp multi-character author chain IDs to single-character PDB chain IDs.
- Preserve the actual residue numbering and CONECT records only when using a
  mature parser; do not text-munge mmCIF.
- Use the original structural document whenever possible. Chain extraction from
  a large complex is acceptable when it is the selected template for the target
  protein, but record the internal chain ID used.

### AlphaFold DB retrieval notes

- Use the AFDB API endpoint:
  `https://alphafold.ebi.ac.uk/api/prediction/{UniProtAccession}`
- The response contains `globalMetricValue` (mean pLDDT, 0-100) and a `pdbUrl`.
- Apply the AFDB quality gate using both available evidence: reject a model if
  the API `globalMetricValue` is below 60, and verify the downloaded PDB's mean
  CA B-factor pLDDT. Direct-download logs must store the actual value used.
- Quarantine or explicitly report pre-existing AFDB-assigned final structures
  that fail the same gate; the wrapper does this by default after direct AFDB
  materialization.
- If a protein maps to multiple UniProt IDs, evaluate all candidates and keep
  one deterministic best model by sequence/mapping evidence.
- If the first UniProt candidate has no AFDB entry, try the next candidate.
- Save the final model as `structure/<protein_sequence_id>.pdb`.
- AFDB files may be downloaded directly as PDB or as mmCIF depending on the
  API response. Verify the downloaded document before conversion or use.
- Treat a `NO_PASSING_AFDB` result as final for that accession and version, but
  make the rule explicit in the log: no model existed, or a mapped model did not
  correspond to the original sequence at an acceptable identity/coverage level.

### Homology Modeling

- Use MODELLER only when the user's pipeline accepts it, the dependency/license
  is available, and the target qualifies for homology modeling.
- Respect all applicable legal and license terms when using MODELLER.
- Build from the exact selected template chain, not the whole PDB complex.
- By default, model one deterministic best qualifying template per protein. If
  the user explicitly asks for 'model every candidate', preserve every
  qualifying candidate and attempt each independently, recording the outcome
  separately before final selection.
- Record template ID, template chain, identity, query coverage, subject
  coverage, alignment boundaries, model DOPE/GA341 or other available metrics,
  and the final output path.
- MODELLER is normally single-threaded per model job. Accelerate large manifests
  by running multiple independent MODELLER processes under the user's bounded
  worker limit; multi-threading inside a single MODELLER run usually provides no
  speedup.
- Preserve separate manifests for PDB-only and AFDB-only homology candidates.
  Use a manifest-specific log prefix or run identifier when they may run
  concurrently, and merge only final provenance afterward.
- Do not run direct structure download and homology modeling concurrently for
  the same target protein set. They may both write `structure/<protein_id>.pdb`.
  Run direct PDB/AFDB materialization first, then homology modeling; or use an
  explicit sequencing barrier before starting any modeling job.
- When writing a final structure file, always check again immediately before
  move/rename that the target does not already contain another non-empty
  structure. If another writer created it, do not overwrite; record `EXISTS` or
  `SKIP_RACE` and keep the already-existing structure.
- Keep the URL actually used for each download in provenance. If a downloaded
  payload was gzip-compressed, record that decompression happened.
- For AFDB fallback/API models, a valid accession can return HTTP 404 if no
  model is currently available. Classify it as `AFDB_404`, not as a transient
  network failure. Recompute or repair the candidate only when a new local AFDB
  database snapshot or mapping table is supplied.
- If homology modeling is unavailable, keep the query `PENDING` and route it to
  Phase B only if that matches the user's requested fallback.

### Phase B: Structure Prediction

Purpose:

- Predict structures only for protein sequences that did not obtain a known
  structure in Phase A.

Rules:

- Build a Phase B pending list by comparing existing
  `structure/<protein_sequence_id>.pdb` files against all input sequence IDs.
- Restrict this list to proteins with no passing PDB, AFDB, or homology-model
  result in the previous steps.
- Predict only proteins whose final structure file does not exist yet.
- Filter predicted structures by `pLDDT >= 0.6` (equivalently mean confidence
  >= 60 if the API reports 0-100). Require at least 60% of residues to pass
  the same 60 threshold, not only a passing mean.
- Save the best passing model as `structure/<protein_sequence_id>.pdb`.
- Log successes, low-confidence results, API errors, and missing results so the
  run can be resumed.

Predictor selection:

- **Do not assume a predictor.**
- Ask the user or use the user-specified option: NVIDIA NIM OpenFold2,
  ESMFold, or another available predictor.
- If the user supplies an NVIDIA NIM OpenFold2 API key or endpoint, use NVIDIA
  NIM OpenFold2. If the user supplies an ESMFold CPU deployment, use that.
- Predictor choice is a user-provided execution policy. Do not silently switch
  predictors or install a large framework without confirmation.
- Record the selected predictor, endpoint or execution environment, model
  version, and whether templates/MSAs were used.

### OpenFold2 API notes

- Endpoint example:
  `https://health.api.nvidia.com/v1/biology/openfold/openfold2/predict-structure-from-msa-and-template`
- Minimal request can use only:
  ```json
  {
    "sequence": "<amino acid sequence>",
    "selected_models": [1, 2],
    "relax_prediction": false
  }
  ```
- The response JSON contains `structures_in_ranked_order`, where each entry
  includes `structure` (PDB text), `confidence` (0-100 mean score), `plddt`
  (per-residue list), and `model_param_set`.
- Keep the highest-confidence structure that passes `confidence >= 60`,
  and keep only structures where at least 60% of residues also have per-residue
  confidence >= 60.
- Store API keys outside the workspace or in a protected,
  user-authorized credentials file; do not commit them.

## Coverage Summary

For each input FASTA file:

- Total proteins.
- Proteins with a PDB experimental structure.
- Proteins with an AlphaFold DB structure.
- AFDB structures additionally split into passing/failing the `pLDDT >= 60`
  gate in the AlphaFold DB retrieval notes.
- Proteins with a homology model.
- Proteins with a Phase B predicted structure.
- Proteins still pending because of dependencies, user policy, or execution
  limits.
- Proteins still missing a structure.

Also report species-level totals using the default input rule above.

Add a workflow-level summary:

- Qualifying FASTA files and excluded/non-empty exceptions.
- Total original query IDs.
- PDB direct-reuse, AFDB direct-retrieval, homology-model, Phase B prediction,
  and missing counts.
- Pending counts due to unavailable predictor/deps or server limits.
- Per-source FASTA counts as machine-readable TSV.

The coverage report must be regenerateable from provenance logs and existing
structure files without rerunning retrieval.
