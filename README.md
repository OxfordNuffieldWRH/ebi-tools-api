# ebi-tools-api
Utility module to query the EBI Tools REST API

Only BLASTp endpoint is implemented, but the codebase may be useful to add other endpoints easily.

### Install

```
pip install git+https://github.com/OxfordNuffieldWRH/ebi-tools-api
```

### Usage


```python
from ebi_tools_api import EBITools

ebi = EBITools(email='your@email.com')

mouse_taxonomy_id = 10090
human_il6 = 'MNSFSTSAFGPVAFSLGLLLVLPAAFPAPVPPGEDSKDVAAPHRQPLTSSERIDKQIRYILDGISALRKETCNKSNMCESSKEALAENNLNLPKMAEKDGCFQSGFNEETCLVKIITGLLEFEVYLEYLQNRFESSEEQARAVQMSTKVLIQFLQKKAKNLDAITTPDPTTNASLLTKLQAQNQWLQDMTTHLILRSFKEFLQSSLRALRQM'

result = ebi.blastp(
    sequence=human_il6,
    title='Search for human IL6 homologues in mouse',
    taxids=mouse_taxonomy_id
)
# return a DataFrame
result.hits
```

Simplified summary:

```python
result.hits_simple
```

Underlying JSON data converted to Python dict:

```python
result.data
```

Nice graphical summary for Jupyter notebooks:

```python
result.summarize()
```

### Caching

By default the results will be cached on disk in `.ebi_tools_cache` to minimise traffic to EBI tools.
To clear cache, remove the contents of this hidden directory. You can change the path to cache in `EBITools` constructor:

```python
ebi = EBITools(email=my_email, cache_dir='some/path/to/cache/dir')
```

### Timeout

After sending an initial query, the tool will check for the result every `i` seconds, where `i` increases by one after each failed attempt.
By default 100 attempts will be made. This can be customised in the constructor:


```python
ebi = EBITools(email=my_email, attempts_threshold=20)
```
