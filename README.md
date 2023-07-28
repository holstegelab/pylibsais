A Python C module wrapper for the [libsais algorithm implementation by Ilya Grebnov](https://github.com/IlyaGrebnov/libsais).
Based on the [PySAIS module by Alexey Gritsenko](https://github.com/AlexeyG/PySAIS).

Extended with functions to speed-up discovery and selection of repeating and/or reoccurring k-mers.



Installation:
---------
```
./setup.py build
./setup.py install

```


Functions:
------------
- sais: construct suffix and LCP array
- min_string: returns the least lexicographic cyclical rotation index of a string (Booth's algorithm)
- max_suffix: returns the maximum suffix length at each string position, bounded by the strings end, or the symbol $ and #. 
- kmer_count: constructs a list of all (non-repetitive) k-mers and their counts from the suffix and LCP array.
- kmer_mask_potential: determines for a k-mer the total masked sequence, max masked per string segment bounded by $, and max continuous masked segement.
- kmer_mask: masks a string with k-mer, potentially limiting to consecutive regions. Returns masked string + string positions that were masked.
- kmer_mask_simple: masks a string with k-mer without using suffix array info. Returns masked string + string positions that were masked. 
- kmer_repetitive: determines if a kmer consists of repetitive sequences. 


Example:
------------
See example.py

```
python example.py example/example_sequence.fa
```

License:
---------
The code is available under the MIT license (see LICENSE file).

