# EvoTrees

Plant phylogenetic trees based on specific genes

# Requirements

- [muscle](https://anaconda.org/bioconda/muscle)
- [biopython](https://biopython.org/)
- [RAxML](https://github.com/stamatak/standard-RAxML)

# Run

```bash
./phy.py demo/species.txt
```

Trees, fastas and other data will be under directory `results` after execution.

# Compare trees

## Requirements

To compare trees, module `robinson-foulds` is required:
- [robinson-foulds](https://pypi.org/project/robinson-foulds/)

But after installation, script should be edited as follow:

```python
t1 = Tree(args.treefile1, format=1)
t2 = Tree(args.treefile2, format=1)
  ```

## Run

``rf 1/RAxML_bestTree.results 2/RAxML_bestTree.results``

Robinson-Foulds distance should be displayed as single float.

# Example

```bash
# plant trees
./phy.py demo/species.txt

# robinson-foulds distance
rf demo/species.nwk results/RAxML_bestTree.results
0.0
```

<html>
<body>
    <div>
        <h4>RAxML_bestTree.results</h4>
        <p>
            <img src="demo/species_raxml_tree.png">
        </p>
        <h4>Real tree from <a href="http://timetree.org/">TIMETREE</a></h4>
        <p>
            <img src="demo/species_tree.png">
        </p>
    </div>
</body>
</html>

