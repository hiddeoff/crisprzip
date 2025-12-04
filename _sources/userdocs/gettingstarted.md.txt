# Getting started


## Installation
CRISPRzip is available on [PyPi](https://pypi.org/) and can be installed 
with [pip](https://pip.pypa.io/en/stable/).

### User installation
Although CRISPRzip can be directly installed from pip, creating a virtual 
environment makes it easier to manage dependencies. Below you can see some 
intructions on how to generate a virtual environment with 
[venv](https://docs.python.org/3/library/venv.html) assuming you already 
have python installed in a bash-like terminal.

```shell
python -m venv crisprzip-venv
source crisprzip-venv/bin/activate
pip install crisprzip
```

### Developer installation
To be able to make changes and contributions to CRISPRzip, you will need to 
get your own copy of the source code and install software dependencies on 
your own. Assuming you have a python and git installed in a bash-like 
terminal, the installation process can be done with the following 
instructions.

```shell
git clone https://github.com/hiddeoff/crisprzip.git
cd crisprzip
python -m venv crisprzip-venv
source crisprzip-venv/bin/activate
pip install -e .
```

Please use our 
[contributing guidelines](https://github.com/hiddeoff/crisprzip/blob/main/CONTRIBUTING.md)
 if you would like us to consider your developments in a future CRISPRzip 
release

### Verifying installation
CRISPRzip development includes a cross-platform compatibility and test 
[workflow](.github/workflows/compatibility-test.yml). If you would like to 
verify your local installation, please follow the Developer Installation 
instructions, then run the following test.

```shell
source crisprzip-venv/bin/activate
pip install -e '.[tests]'  # installs pytest and pandas
pytest tests/cleavage_binding_prediction/test_cleavage_binding_prediction.py -v
```

You can also follow the user installation and execute the code in the 
Usage section below.

## Usage
CRISPRzip makes predictions about cleavage and binding activity on on- and
off-targets. First, you define the protospacer and target sequence, and then,
you can predict the fraction cleaved or bound.

```python
# 1. load parameter set
from crisprzip.kinetics import load_landscape
searcher = load_landscape("sequence_params")

# 2. define Cas9, gRNA and DNA target
searchertargetcomplex = searcher.probe_sequence(
    protospacer = "AGACGCATAAAGATGAGACGCTGG",
    target_seq  = "AGACCCATTAAGATGAGACGCGGG",  # A13T G17C
)

# 3. predict activity
f_clv = searchertargetcomplex.get_cleaved_fraction(
    time=600,  # 10 minutes
    on_rate=1E-1
)
f_bnd = searchertargetcomplex.get_bound_fraction(
    time=600,  # 10 minutes
    on_rate=1E-1
)

# 4. format output
print(f"After 10 minutes, the target (A13T G17C) is ...")
print(f"- cleaved for {100 * f_clv:.1f}% by Cas9")
print(f"    or  ")
print(f"- bound for {100 * f_bnd:.1f}% by dCas9")
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
