# mad4hatter_utilities
Compilation of scripts to work with mad4hatter's outputs


## Filtering 💩🧹



## Format resmarkers 💅💇‍♀️✨

`format_resmarkers.R` converts mad4hatter's resmarker_table.txt (long format) into the former "report" format. Read counts are shown between brackets. Note: when there are more than one resistance marker for a given position, the reference variant (REF) is displayed first.

### Usage

```shell
Rscript format_resmarkers.R --input resmarkers_table.txt --output resmarkers_table_old_format.csv
```

## Analyses and visualizations 📊 👀
