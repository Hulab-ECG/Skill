# Known-Structure Retrieval

This reference expands Stages 2-5. Use only the parts matching the selected
execution route. Do not repeat both AFDB routes for the same protein unless the
user asks for a re-selection pass.

## Route Selection

1. **Local sequence search is the default route.** It searches the user-provided
   local PDB sequence library and a local AFDB full-sequence FASTA library
   exactly like other template-search workflows.
2. **AFDB API retrieval is a fallback**, used when there is no local AFDB
   sequence library or when the user explicitly asks for a mapping-based
   supplement.
3. A protein removed from pending by a passing PDB or AFDB result must not be
   sent to another known-structure route.

Record the selected route, database version/snapshot, FASTA path, DIAMOND
version, checksums, sequence counts, and user-specified thread/run limits.

## Reference Databases

### PDB

- Accept an explicit `--pdb-db-dir` or equivalent configuration.
- A ready-to-use library can be `*.dmnd` or `pdb_seqres.txt` /
  `pdb_seqres.txt.gz`.
- If only FASTA exists, create the DIAMOND database in the same user-specified
  directory with `diamond makedb`.
- Never invent or overwrite a database directory. Ask for a path or obtain
  explicit permission before creating one.
- Validate that required files are non-empty, PDB headers are present, and the
  DIAMOND database is readable by the installed DIAMOND version.
- Rebuild only when the library is missing, unreadable, stale, or the user asks
  for a refresh. Use a filesystem lock for long builds.

### AFDB

- An explicit AFDB FASTA file, such as a full-sequence archive, is built into a
  DIAMOND database like the PDB library.
- Verify download completeness before indexing; do not rely on file size
  alone. Validate the FASTA terminator and record the sequence count.
- Store the AFDB version/snapshot and database checksum. This is especially
  important because AFDB models are versioned and updated over time.

## Sequence Search

Canonical search parameters:

```bash
diamond blastp --ultra-sensitive --query-cover 50 --subject-cover 50 --id 50 --iterate
```

- Run searches separately against PDB and AFDB so hits remain distinguishable.
- Save all threshold-passing hits, not just the first/best hit.
- Preserve `pident`, `qcovhsp`, `scovhsp`, e-value/score, subject ID, and the
  original query ID.
- DIAMOND only performs candidate discovery. Final candidate selection is done
  separately in Stages 5 and 8.

Hit classes:

1. **Direct reuse**: 100% sequence identity and full bidirectional coverage.
   This is the only local-search class eligible for direct structure reuse.
2. **Template-only**: passes the search thresholds but is not exact or not
   full length. Do not reuse directly; route by Stage 5.

"Full bidirectional coverage" must be represented explicitly, for example as
both `qcovhsp >= 100` and `scovhsp >= 100`, or by an equivalent exact
whole-chain alignment definition if the DIAMOND columns are not both used.

## PDB Structural Materialization

For every direct-reuse PDB chain:

1. Parse the subject ID. Example: `1abc_A` gives `pdb_id=1abc`, `chain=A`.
2. Download `https://files.rcsb.org/download/{PDB_ID}.pdb` first.
3. If legacy `.pdb` returns 404, retry
   `https://files.rcsb.org/download/{PDB_ID}.cif`.
4. Parse with Bio.PDB, gemmi, or another mature parser. Never text-munge mmCIF.
5. Extract the selected chain. Do not copy the whole complex unless the matched
   target is the whole entry.
6. Handle multi-character author chain IDs with the parser.
7. Write the final file atomically to
   `<species_dir>/structure/<original_protein_sequence_id>.pdb`.
8. Record PDB ID, chain, method, resolution when available, download status,
   provenance, and final path.

PDB preference ordering is X-ray/crystallography > Cryo-EM > NMR > other.
Within the same method, prefer lower resolution. Keep only one best PDB model
per original protein ID.

mmCIF may remain as a cache but the final per-protein structure must use the
`.pdb` filename required by the project.

## AFDB Local Search Route

For an exact AFDB hit:

1. Resolve the subject ID to the AFDB model accession.
2. Always fetch the official API response once to obtain the exact `entryId`
   and `pdbUrl` before downloading. Do not guess static `model_v1`, `model_v4`,
   or `model_v6` URLs.
3. Verify the downloaded document is complete and parseable. If the `pdbUrl`
   payload is gzip-compressed, decompress it before validation.
4. Save one passing model per original protein ID.

For non-exact AFDB hits, pass them to the template-evaluation stage. They may
be eligible only for homology modeling; they are not direct structure reuse.

## AFDB API Route

Endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/{UniProtAccession}
```

- The response contains `globalMetricValue` (mean pLDDT, 0-100) and model
  download URLs.
- Evaluate all mapped UniProt candidates, not only the first.
- Keep the best passing model and save it under the original protein ID.
- If the first candidate returns 404, try the next candidate.
- A 404 is a valid "AFDB has no model" result for that accession/version,
  provided the request was well formed. Treat 4xx transport/auth failures,
  429, 5xx, timeouts, malformed responses, and incomplete downloads as retry-
  able rather than as valid absence.
- A `NO_PASSING_AFDB` result must explicitly say whether no model existed, all
  models were below the applicable sequence criterion, or the mapped model did not correspond to the
  query at an acceptable identity/coverage level.

AFDB can complement, but must not override, a valid PDB structure.

## UniProt Mapping

### Primary mapping

- RefSeq-style `WP_`, `NP_`, and `YP_` accessions usually work through the
  `RefSeq_Protein` idmapping channel.
- INSDC-style IDs such as `CAP`, `MFC`, and `CAQ` are often rejected by that
  channel.

### Recovery

Use genome annotation or sequence identity bridges:

- Extract xrefs from `gbff`/`gff3`, including `RefSeq:WP_...`,
  `UniParc:UPI...`, `UniRef:UniRef100_...`, `UniRef:UniRef50_...`, and
  `UniRef:UniRef90_...`.
- Use UniParc or sequence checksums to recover possible UniProt candidates.
- Store all candidate UniProt accessions and their provenance. Do not discard
  one-to-many ambiguity; the best passing AFDB model may come from any
  candidate.
- Merge UniParc-derived candidates back into the canonical protein-to-UniProt
  mapping table.

### Batch API notes

- UniProt idmapping stream jobs do not support an extra `fields` parameter in
  the batch submit request; adding `?fields=` can return an empty
  `{results: []}`.
- The old `/idmapping/status/{jobId}` endpoint may return 404; use the current
  job status/result polling workflow and allow long waits for large batches.
- Split very large inputs into chunks, write mapped chunks separately, merge
  later, and retain chunk-level retry state.
