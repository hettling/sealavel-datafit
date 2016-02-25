# !/bin/bash
(printf 'compaction\t' \
	&& (echo `head -1 sealevel-comp-min.tsv` | tr " " "\t") \
	&& printf "min\t" && tail -1 sealevel-comp-min.tsv \
	&& printf "avg\t" && tail -1 sealevel-comp-avg.tsv \
	&& printf "max\t" && tail -1 sealevel-comp-max.tsv) \
> results.tsv
