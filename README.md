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
| sample_1  | V [1968]   | C/V/M/N/K [2872]   | N/I [2151]   | A/I/R [1950]   | T [564]  | I/S/G [3205]     |
| sample_2    | V [838]    | C/V/M/N/K [1245]   | N/I [997]    | A/I/R [845]    | T [127]  | I/S/G [1499], I/S/A [161] |
| sample_3  | V [2516]   | C/V/M/N/K [2754]   | N/I [2426]   | A/I/R [2191]   | T [911]  | I/S/G [2878]     |
| sample_4  | V [1121]   | C/V/M/N/K [1130]   | N/I [1126]   | A/I/R [1008]   | T [328]  | I/S/G [1379]     |
| sample_5  | V [1601]   | C/V/M/N/K [2091]   | N/I [1605]   | A/I/R [1403]   | T [460]  | I/S/G [373], I/A/A [1753] |
| sample_6  | V [1316]   | C/V/M/N/K [1633]   | N/I [1261]   | A/I/R [1425]   | T [479]  | I/S/G [1741]     |
| sample_7  | V [2566]   | C/V/M/N/K [2920]   | N/I [2246]   | A/I/R [2293], A/N/R [137] | T [758]  | I/S/G [2811]     |
| sample_8  | V [2920]   | C/V/M/N/K [2809]   | N/I [2760]   | A/I/R [1942], A/N/R [348] | T [1147] | I/S/G [2784]     |
| sample_9  | V [2626]   | C/V/M/N/K [2571]   | N/I [2607]   | A/I/R [2007]   | T [1462] | I/S/G [2608]     |
| sample_10  | V [701]    | C/V/M/N/K [767]    | N/I [566]    | A/I/R [649]    | T [298]  | I/S/G [651], I/S/A [125] |
| sample_11  | V [1406]   | C/V/M/N/K [1801]   | N/I [1436]   | A/I/R [1223]   | [NA]     | I/S/G [1814]     |



## Analyses and visualizations 📊 👀
