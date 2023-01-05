from typing import Literal, Union, Callable, Coroutine
from dataclasses import dataclass
from functools import cached_property, partial
from json import loads
from pathlib import Path
import asyncio
import pickle
import hashlib

from pandas import DataFrame, isnull
from IPython.display import SVG, Image, display
import requests


def one(values, context):
    if len(values) > 1:
        print(f'Multiple values found for {context}')
    return values[0]


graphics_wrappers = {
    'svg': SVG,
    'png': partial(Image, format='png'),
    'jpg': partial(Image, format='jpg')
}

protein_evidence_of_existence = {
    1: '1. Experimental evidence at protein level',
    2: '2. Experimental evidence at transcript level',
    3: '3. Protein inferred from homology',
    4: '4. Protein predicted',
    5: '5. Protein uncertain'
}

from contextlib import contextmanager

class CacheResource:
    def __init__(self):
       self.result = None

    def set_result(self, result):
       self.result = result


@contextmanager
def disk_cache(self, endpoint: str, cached_only=False, **kwargs):
    reset_cache = kwargs.pop('reset_cache', False)
    verbose = kwargs.pop('verbose', False)
    kwargs = {
      k: kwargs[k]
      for k in sorted(kwargs)
    }
    params_hash = hashlib.sha512(str(kwargs).encode()).hexdigest()
    label = params_hash  # TODO: find a better way?
    cache_path: Path = self.cache_dir / endpoint / str(params_hash)
    if verbose and reset_cache:
        print(
            f'Reseting cache for {label} as requested'
            if cache_path.exists() else
            f'Cache reset for {label} requested but there was no cache'
        )
    resource = CacheResource()
    if cache_path.exists() and not reset_cache:
        if verbose:
            print(f'Cached result found for {label}, loading from disk')
        with open(cache_path, 'rb') as f:
            data = pickle.load(f)
        # guard against collisions
        assert data['kwargs'] == kwargs
        resource.set_result(data['result'])
        yield resource
    else:
        if cached_only:
            raise Exception(f'Cached copy not found in {cache_path}')
        if verbose:
            print(f'{cache_path} not found fetching...')
        yield resource
        data = resource.result
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, 'wb') as f:
            pickle.dump({
                'result': data,
                'kwargs': kwargs
            }, f)


def disk_cached(func: Callable[..., Coroutine]):
    async def wrapper(self, endpoint: str, cached_only=False, **kwargs):
        with disk_cache(self, endpoint, cached_only, **kwargs) as cache_resource:
            if cache_resource.result:
                return cache_resource.result
            else:
                data = await func(self, endpoint, **kwargs)
                cache_resource.set_result(data)
                return data
    return wrapper


def disk_cached_sync(func: Callable):
    def wrapper(self, endpoint: str, cached_only=False, **kwargs):
        with disk_cache(self, endpoint, cached_only, **kwargs) as cache_resource:
            if cache_resource.result:
                return cache_resource.result
            else:
                data = func(self, endpoint, **kwargs)
                cache_resource.set_result(data)
                return data
    return wrapper


@dataclass
class Result:
    job_id: str
    tool: str
    server: str = 'https://www.ebi.ac.uk/Tools/services/rest'
    cache_dir: Path = Path('.ebi_tools_cache')

    def result_types(self, **kwargs):
        return requests.get(
            f'{self.server}/{self.tool}/resulttypes/{self.job_id}',
            **kwargs
        )

    def _get(self, output: str, **kwargs):
        return self._get_cached(
            endpoint=self.tool,
            output=output,
            job_id=self.job_id,
            **kwargs
        )

    @disk_cached_sync
    def _get_cached(self, endpoint: str, output: str, job_id: str, **kwargs):
        return (
            requests.get(
                f'{self.server}/{endpoint}/result/{job_id}/{output}',
                **kwargs
            )
        )


@dataclass
class JSONResult(Result):

    @property
    def program(self):
        return self.data['program']

    @property
    def version(self):
        return self.data['version']

    @cached_property
    def data(self) -> dict:
        return loads(self._get('json').text)


def _int_or_none(x):
    if isnull(x):
        return None
    return int(x)


@dataclass
class BlastResult(JSONResult):
    job_id: str
    tool: str = 'ncbiblast'
    image_format: Literal['svg', 'png', 'jpg'] = 'svg'

    @property
    def hits(self) -> DataFrame:
        df = DataFrame(self.data['hits'])
        if df.empty:
          return df
        return df.set_index('hit_num')

    @staticmethod
    def simplify(hits: DataFrame) -> DataFrame:
        return (
            hits
            .assign(
                existence=lambda df: df.hit_uni_pe.map(_int_or_none).map(protein_evidence_of_existence),
                database=lambda df: df.hit_db.map({
                    'SP': 'Swiss-Prot (Reviewed)',
                    'TR': 'TrEMBL (Unreviewed)'
                }),
                uniprot_link=lambda df: ('<a href="' + df.hit_url + '" target="_blank">' + df.hit_acc + '</a>')
            )
            .drop(columns=[
                'hit_def', 'hit_hsps', 'hit_desc',
                'hit_xref_url', 'hit_url',
                'hit_uni_ox', 'hit_db',
                'hit_uni_pe'
            ])
            .rename(columns=lambda column: column.replace('hit_', ''))
            .rename(columns={
                'uni_de': 'description',
                'uni_os': 'species',
                'uni_gn': 'gene_name',
                'uni_sv': 'sequence_version',
                'len': 'length',
                'acc': 'accession',
                'id': 'identifier'
            })
        )

    @staticmethod
    def extract_local_alignment(
        hits: DataFrame,
        sort: str = 'hsp_expect',
        ascending: bool = False
    ):
        return hits.apply(
            lambda row: (
                DataFrame(row.hit_hsps)
                .sort_values(sort, ascending)
                .iloc[0]
            ),
            axis=1
        ).rename(columns=lambda column: column.replace('hsp_', ''))

    def hits_simple(self) -> DataFrame:
      return self.simplify(self.hits)

    @cached_property
    def visual(self) -> Union[SVG, Image]:
        return self._wrap_image(self._get(f'visual-{self.image_format}').content)

    @cached_property
    def fast_family_and_domain_prediction(self) -> Union[SVG, Image]:
        return self._wrap_image(self._get(f'ffdp-subject-{self.image_format}').content)

    def _wrap_image(self, data):
        image_wrapper = graphics_wrappers[self.image_format]
        return image_wrapper(data=data)

    @cached_property
    def _by_accession(self):
        return self.hits_simple().set_index('accession')

    def __getitem__(self, accession):
        return {
            'chosen': {
                'accession': accession,
                **self._by_accession.loc[accession].to_dict(),
            },
            'all_results': self
        }

    @property
    def hits_summary(self):
        return self.hits_simple().style.format({'e_value': '{:e}', 'identity': '{:.2f}'})

    def summarize(self):
        display(self.hits_summary)
        display(self.visual)
        display(self.fast_family_and_domain_prediction)

    def __repr__(self):
        return f'<{self.version} with {len(self.hits)} results>'


def parse_needle_summary(text: str):
    data = {}
    for line in text.split('\n'):
        if line.startswith('#') and ':' in line:
            key, value = [x.strip() for x in line[1:].split(':', maxsplit=1)]
            data[key] = value
    return data


@dataclass
class NeedleResult(Result):
    job_id: str
    tool: str = 'emboss_needle'

    @property
    def out(self) -> dict:
        text = self._get('out').text
        return parse_needle_summary(text)

    @property
    def alignment(self) -> dict:
        return self._get('aln').text


@dataclass
class StretcherResult(Result):
    job_id: str
    tool: str = 'emboss_stretcher'


EBI_HEADERS = {
    'Accept': 'text/plain',
    'Content-Type': 'application/x-www-form-urlencoded'
}


@dataclass
class EBITools:
    email: str
    server: str = 'https://www.ebi.ac.uk/Tools/services/rest'
    cache_dir: Path = Path('.ebi_tools_cache')
    attempts_threshold: int = 50
    backoff_limit: int = 10
    verbose: bool = True

    async def blastp(self, sequence, exp='1e-10', stype='protein', database='uniprotkb', **query) -> BlastResult:
        job_id = await self._query_cached(
            endpoint='ncbiblast',
            program='blastp',
            task='blastp',
            exp=exp,
            stype=stype,
            database=database,
            sequence=sequence,
            **query
        )
        return BlastResult(
            job_id=job_id,
            server=self.server,
            cache_dir=self.cache_dir
        )

    async def needle(self, asequence, bsequence, matrix='EBLOSUM62', stype='protein', database='uniprotkb', **query) -> NeedleResult:
        """EMBOSS Needle creates an optimal global sequence alignment of two input sequences using the Needleman-Wunsch alignment algorithm."""
        if not asequence:
            raise ValueError('First sequence is missing')
        if not bsequence:
            raise ValueError('Second sequence is missing')
        job_id = await self._query_cached(
            endpoint='emboss_needle',
            matrix=matrix,
            stype=stype,
            database=database,
            asequence=asequence,
            bsequence=bsequence,
            **query
        )
        return NeedleResult(
            job_id=job_id,
            server=self.server,
            cache_dir=self.cache_dir
        )

    async def stretcher(self, asequence, bsequence, matrix='EBLOSUM62', stype='protein', database='uniprotkb', **query) -> StretcherResult:
        """EMBOSS Stretcher creates an optimal global sequence alignment of two input sequences using O(min(N, M)) space."""
        if not asequence:
            raise ValueError('First sequence is missing')
        if not bsequence:
            raise ValueError('Second sequence is missing')
        job_id = await self._query_cached(
            endpoint='emboss_stretcher',
            matrix=matrix,
            stype=stype,
            database=database,
            asequence=asequence,
            bsequence=bsequence,
            **query
        )
        return StretcherResult(
            job_id=job_id,
            server=self.server,
            cache_dir=self.cache_dir
        )

    @disk_cached
    async def _query_cached(self, endpoint: str, **query) -> str:
        params = {'email': self.email, **query}
        query = requests.post(
            f'{self.server}/{endpoint}/run',
            params=params,
            headers=EBI_HEADERS
        )
        if query.status_code != 200:
            raise ValueError(f'Status code: {query.status_code}, content: {query.content}')
        job_id = query.text

        status = None
        i = 1

        while status != 'FINISHED':
            if i > self.attempts_threshold:
                raise Exception(f'Not finished at {self.attempts_threshold}th attempt')
            await asyncio.sleep(i)
            status = (
                requests.get(
                    f'{self.server}/{endpoint}/status/{job_id}',
                    headers={'Accept': 'text/plain'},
                )
                .text
            )
            assert status in {'RUNNING', 'FINISHED'}, status
            if self.verbose:
                print('.', end='')
            if i < self.backoff_limit:
                i += 1

        return job_id

