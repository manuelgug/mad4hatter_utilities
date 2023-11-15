# mad4hatter_utilities

Compilation of scripts to work with mad4hatter's outputs

## Filtering 💩🧹

`filtering.R` applies contaminants filtering and minimum allele frequency (MAF) filtering to remove potential false positives from allele_data.txt, resmarker_microhap_table.txt and resmarkers_table.txt based on read counts from negative controls. It also appends allele counts and frequencies to the aforementioned inputs.

### Usage

```bash
Rscript filtering.R [--allele_table PATH] [--microhaps_table PATH] [--resmarker_table PATH] [--CFilteringMethod METHOD] [--MAF VALUE] [--exclude_file PATH] [--use_case_amps PATH] [--outdir PATH]
```

- `--allele_table`: Path to the input allele table.
  
- `--microhaps_table`: Path to the input microhaps table (optional, but required if `--resmarkers_table` is provided).

- `--resmarkers_table`: Path to the input resmarkers table (optional).

- `--CFilteringMethod`: Contaminants filtering method. Options: *global_max*, *global_q95*, *amp_max*, *amp_q95*. Default: *global_max*.

  - *global_max*: Single threshold derived from the maximum read count from all amplicons across negative controls.
  - *global_q95*: Single threshold derived from the 95th percentile of read counts from all amplicons across negative controls.
  - *amp_max*: Amplicon-specific thresholds derived from the maximum read count for each amplicon across negative controls.
  - *amp_q95*: Amplicon-specific thresholds derived from the 95th percentile of read counts for each amplicon across negative controls.

- `--MAF`: Minimum allele frequency filter. Default: 0.

- `--exclude_file`: Path to the file containing sampleIDs to exclude (optional).
  
- `--use_case_amps`: Path to the file containing IDs of the use case amplicons (optional).

- `--outdir`: Name of the directory where results will be stored.

### Nomenclature of controls

- **Positive Controls**: The script identifies positive controls using the following conditions:
  - Sample IDs containing "3D7" (case-insensitive) and not containing "Dd2", "HB3", or "PM" (case-insensitive).

- **Negative Controls**: The script identifies negative controls using the following conditions:
  - Sample IDs containing "Neg" (case-insensitive).

## Format resmarkers 💅💇‍♀️✨

`format_resmarkers.R` converts mad4hatter's resmarker_table.txt and/or resmarker_microhap_table.txt (long format) into their former "report" format. 

Read counts are shown between brackets and multiple variants are separated by commas.

Note: when there are multiple resistance markers or haplotypes for a given position, the reference variant (REF) is displayed first.

### Usage

```shell
Rscript format_resmarkers.R --input_res resmarkers_table.txt --output_res resmarkers_table_old_format.csv --input_haps resmarkers_microhap_table.txt --output_haps resmarkers_microhap_table_old_format.csv
```
### Example output (haplotypes)

|              | arps10_127 | crt_72/73/74/75/76 | dhfr_108/164 | dhfr_16/51/59 | dhfr_185 | dhps_431/436/437 |
|--------------|------------|---------------------|--------------|--------------|----------|------------------|
| S1  | V [1601]   | C/V/M/N/K [2091]   | N/I [1605]   | A/I/R [1403]   | T [460]  | I/S/G [373], I/A/A [1753] |
| S2  | V [1316]   | C/V/M/N/K [1633]   | N/I [1261]   | A/I/R [1425]   | T [479]  | I/S/G [1741]     |
| S3  | V [2566]   | C/V/M/N/K [2920]   | N/I [2246]   | A/I/R [2293], A/N/R [137] | T [758]  | I/S/G [2811]     |
| S4  | V [701]    | C/V/M/N/K [767]    | N/I [566]    | A/I/R [649]    | T [298]  | I/S/G [651], I/S/A [125] |
| S5| V [1406]   | C/V/M/N/K [1801]   | N/I [1436]   | A/I/R [1223]   | [NA]     | I/S/G [1814]     |



## Analyses and visualizations 📊 👀

## K13 non-synonymous mutations screening
