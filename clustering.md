# Summary

100 Pfam familes were chosen randomly, merged and stored in `merged.fasta` file.
Then 2 tools were used to cluster them:
- [mmseqs2](https://github.com/soedinglab/MMseqs2)
- [cd-hit](https://github.com/weizhongli/cdhit)

Results were much more similar to real clusters for `mmseqs2`:

```bash
Clustering results:
	[Jaccard similarity] mmseqs2 = 0.160514
	[Jaccard similarity] cd-hit = 0.017746
```

To obtain the same results, run following command (keep in mind that `cd-hit` and `mmseqs2` should be in your `$PATH` and `merged.fasta` file must be in the same directory):

```bash
./pfam.py
```
