# ebi-tools-api
Utility module to asynchronously query the EBI Tools REST API

Only BLASTp and Needle endpoints are implemented, but the codebase may be useful to add other endpoints easily.

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

result = await ebi.blastp(
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
To clear cache, either:
- pass `reset_cache=True` to a query function (e.g. `blastp(sequence=human_il6, reset_cache=True)` **preserving all other arguments as in original call**
- remove the contents of this hidden directory to clear entire cache at once

Pass `verbose=True` to see cache status.

You can change the path to cache in `EBITools` constructor:

```python
ebi = EBITools(email=my_email, cache_dir='some/path/to/cache/dir')
```

### Timeout

After sending an initial query, the tool will check for the result every `i` seconds, where `i` increases by one after each failed attempt,
up to a configurable maximum interval (default 5 seconds). By default 100 attempts will be made. This can be customised in the constructor:


```python
ebi = EBITools(email=my_email, attempts_threshold=20, backoff_limit=10)
```
