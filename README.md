# mad4hatter_utilities
Compilation of scripts to work with mad4hatter's outputs


## Filtering 💩🧹



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
