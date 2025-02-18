# Getting started


## Installation
CRISPRzip is on [PyPi](https://pypi.org/) and can be installed 
with [pip](https://pip.pypa.io/en/stable/).

```shell
pip install crisprzip
```

## Usage
CRISPRzip makes predictions about cleavage and binding activity on on- and
off-targets. First, you define the protospacer and target sequence, and then,
you can predict the fraction cleaved or bound.
```python
# 1. load parameter set
import json
with open('data/landscapes/sequence_params.json', 'r') as file:
    sequence_params = json.load(file)['param_values']

# 2. define Cas9, gRNA and DNA target
from crisprzip.kinetics import *
searchertargetcomplex = SearcherSequenceComplex(
    protospacer = "AGACGCATAAAGATGAGACGCTGG",
    target_seq  = "AGACCCATTAAGATGAGACGCGGG",  # A13T G17C
    **sequence_params
)

# 3. predict activity
f_clv = searchertargetcomplex.get_cleaved_fraction(
    time=600, # 10 minutes
    on_rate=1E-1
)
f_bnd = searchertargetcomplex.get_bound_fraction(
    time=600, # 10 minutes
    on_rate=1E-1
)

# 4. format output
print(f"After 10 minutes, the target (A13T G17C) is ...")
print(f"- cleaved for {100*f_clv:.1f}% by Cas9")
print(f"    or  ")
print(f"- bound for {100*f_bnd:.1f}% by dCas9")
```
Output:
```
After 10 minutes, the target (A13T G17C) is ...
- cleaved for 10.5% by Cas9
    or  
- bound for 94.2% by dCas9
```

See the [tutorial](./tutorial.ipynb) for more examples how to explore 
sequence, time and concentration dependency.
