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

## Analyses and visualizations 📊 👀
