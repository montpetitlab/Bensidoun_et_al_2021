To run:

```
snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=256000 --cluster "sbatch -t 480 -J bensi8 -p med2 -n 1 -N 1 -c 2 --mem=8GB" -k
```
