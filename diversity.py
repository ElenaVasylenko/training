import json
import biom
from skbio.diversity import alpha_diversity, beta_diversity
from os.path import isfile
from biom.table import Table
from biom import load_table
from numpy import asarray

from diversity_data import input_biom



ALPHA_DIVERSITY_METHODS = ['chao1', 'simpson', 'shannon']
BETA_DIVERSITY_METHODS = ['abund_jaccard', 'bray_curtis']


def get_biom_table(data):
    """returns a biom object regardless of whether path or object given"""
    try:
        if isfile(data):
            return load_table(data)
    except TypeError:
        if type(data) == Table:
            return data
        else:
            raise TypeError('Data is neither a path to a biom table or a' +
                            ' biom table object.')


def alpha(table: biom.Table):
    """

    :param table:
    :return:
    """
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    table = get_biom_table(table)
    alpha_diversities = []
    counts = table.matrix_data.toarray().astype(float).T
    sample_ids = table.ids(axis='sample')
    sample_metadata = dict(zip(table.ids(), table.metadata()))

    for metric in ALPHA_DIVERSITY_METHODS:
        result = alpha_diversity(metric=metric, counts=counts, ids=sample_ids)
        result.name = metric
        alpha_diversities.append(result)
    aggregated_diversity_results = aggregate_results(alpha_diversities, sample_ids)
    formatted_diversity_results = _format_alpha_results(aggregated_diversity_results, sample_metadata)

    return formatted_diversity_results


def aggregate_results(table, sample_ids):
    """

    :param table:
    :param sample_ids:
    :return:
    """
    diversity_aggregated = (asarray(table).transpose(), sample_ids, ALPHA_DIVERSITY_METHODS)
    return diversity_aggregated


def _format_alpha_results(result, sample_metadata=None):
    """

    :param result:
    :param sample_metadata:
    :return:
    """
    if sample_metadata is None:
        sample_metadata = {}

    data, sample_names, calc_names = result
    sample_list = []
    labels_mapping = {}

    for sample_index, sample_name in enumerate(sample_names):
        single_sample_metadata = sample_metadata.get(sample_name, {})
        sample_dict = dict(sample_name=sample_name, values=dict(),
                           metadata=dict(label=single_sample_metadata.get('label', 'label0')))
        labels_mapping[single_sample_metadata.get('label', 'label0')] = \
            single_sample_metadata.get('label_name', 'No Label')
        for calc_index, method_name in enumerate(calc_names):
            value = data[sample_index][calc_index]
            if str(value) == 'nan':
                # when sample has 0-hits, alpha diversity returns nan value
                value = 0
            sample_dict['values'][method_name] = value
        sample_list.append(sample_dict)

    return dict(data=sample_list, methods=calc_names, labels=labels_mapping)

analysis_alpha_diversity = alpha(input_biom)
print(analysis_alpha_diversity)
