"""Usage: get_representative.py <pangolin_output> <embl_seqs> [-o OUTFILE]

-o OUTFILE  Output file name [default: repr_seqs.tsv]
"""
import csv
from docopt import docopt
import pandas as pd


def main(args):
    pangolin_output = args.get('<pangolin_output>', 'lineages.csv')
    embl_seqs = args.get('<embl_seqs>', 'seqs_embl_covid18.tsv')

    inclusion_accessions = ['MT929124', 'MZ009837']

    lineages = pd.read_csv(pangolin_output)
    lineages = (
        lineages
            .assign(accession=lineages.taxon.map(extract_accession))
            .set_index('accession')
    )
    require_cols(lineages, ['lineage'])

    seqs = pd.read_table(
        embl_seqs,
        parse_dates=['collection_date'],
        quoting=csv.QUOTE_ALL,
        index_col='accession_id',
    )
    require_cols(seqs, ['collection_date'])

    joined = seqs.join(lineages, how="inner")
    joined = joined.assign(accession_id=joined.index)

    rep_seqs = (joined
        .pipe(create_date_col)
        .pipe(include, inclusion_accessions)
        .pipe(remove_failed)
        .pipe(remove_older_than, 2019_12_20)
        .pipe(exclude_ambiguous_dates)
        .pipe(remove_coverage_under, 98.0)
        .pipe(keep_most_supported_vocs)
        .pipe(earliest_date_per_lineage)
        .pipe(drop_lineage_none)
        .pipe(lambda x: x.sort_values(['scorpio_call', 'lineage']))
    )
    # print(rep_seqs.columns)
    # print(rep_seqs)
    keep_cols = [
        'accession_id',
        'cross_references',
        'date',
        'country',
        'lineage',
        'ambiguity_score',
        'scorpio_call',
        'scorpio_support',
        'scorpio_conflict',
        'coverage',
        'pangolin_version',
        'note',
    ]

    # import pdb; pdb.set_trace()
    rep_seqs[keep_cols].to_csv(args.get('-o'), sep='\t', index=None)


def require_cols(df: pd.DataFrame, required_cols: list):
    valid =  all(col in df.columns for col in required_cols)
    if not valid:
        print(f"Error: Dataframe did not contain "
              f"all required columns: {required_cols}")
        exit(1)


def extract_accession(taxon: str) -> str:
    return taxon.split('|')[1]


def create_date_col(df):
    require_cols(df, ['collection_date'])
    return df.copy().assign(
        date=df.collection_date.map(format_date)
    )


def remove_failed(df: pd.DataFrame):
    require_cols(df, ['status'])
    return df.copy()[df.status.str.contains('passed_qc', na=False)]


def remove_coverage_under(df: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    return df.copy()[df.coverage > cutoff]


def remove_older_than(df: pd.DataFrame, date: int):
    # Remove sequences with a spurious collection date of before the outbreak
    df2 = df.copy()[~df.collection_date.isna()]
    return df2[df2.collection_date.astype(int) > date]


def keep_most_supported_vocs(df):
    voc_rows = ~df.scorpio_support.isna()
    vocs = df.copy()[voc_rows]
    rest = df.copy()[~voc_rows]
    top_vocs = vocs.sort_values('scorpio_support', ascending=False).groupby('lineage').head(1)
    # top_rest = rest.sort_values('ambiguity_score', ascending=False).groupby('lineage').head(1)
    return pd.concat([top_vocs, rest])


def include(df, accessions: list[str]):
    inclusion_rows = df.accession_id.isin(accessions)
    inclusion_df = df.copy()[inclusion_rows]

    lineages_to_remove = df[inclusion_rows].lineage.tolist()
    filtered_df = df.copy()[~df.lineage.isin(lineages_to_remove)]
    return pd.concat([inclusion_df, filtered_df])


def exclude_ambiguous_dates(df):
    return df.copy()[~df.date.str.endswith('XX')]


def earliest_date_per_lineage(df):
    return df.copy().sort_values('date').groupby('lineage').head(1)


def drop_lineage_none(df):
    return df.copy().query("lineage != 'None'")


def format_date(date: str) -> str:
    """
    Accepts an EBI json date string and formats it for Nextstrain
    e.g. 20210100 -> 2021-01-XX
    """
    if date is None: return "?"
    date = str(date)
    year = date[:4]
    month = date[4:6].replace("00", "XX")
    day = date[6:].replace("00", "XX")
    return f"{year}-{month}-{day}"


if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

