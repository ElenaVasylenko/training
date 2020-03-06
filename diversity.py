import json

import biom
import numpy
import skbio
from termcolor import colored
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.diversity.alpha import chao1, shannon, simpson
from os.path import abspath, basename, exists, dirname, join, splitext, isfile
from biom.table import Table
from biom import load_table
from numpy import array


ALPHA_DIVERSITY_METHODS = ['chao1', 'simpson', 'shannon']
BETA_DIVERSITY_METHODS = ['abund_jaccard', 'bray_curtis']


# def alpha_diversity(biom, target):
#     print("Running alpha diversity analysis with methods %s", ALPHA_DIVERSITY_METHODS)
#
#     with open(target, 'w') as file_:
#         data = _single_file_alpha(biom, ALPHA_DIVERSITY_METHODS)
#         file_.write(json.dumps(data))
#

def alpha_diversity(biom, file):
    print("Running alpha diversity analysis with methods %s", ALPHA_DIVERSITY_METHODS)
    data = _single_file_alpha(biom, ALPHA_DIVERSITY_METHODS)

def getBiomData(data):
    """returns a biom object regardless of whether path or object given"""
    try:
        if isfile(data):
            return load_table(data)
    except TypeError:
        if type(data) == Table:
            otu_table = data
            return otu_table
        else:
            raise TypeError('Data is neither a path to a biom table or a' +
                            ' biom table object.')

def _single_file_alpha(biom, metrics):
    """
    Calculate Alpha diversity.

    :param biom: Input BIOM table. Can accept both paths and biom.Table
    :type biom: biom.Table or str
    """
    metrics_list = metrics
    biom_data = getBiomData(biom)
    transposed_biom_data = biom_data.transpose()
    calcs = []

    for metric in metrics_list:
        alpha = []
        c = alpha_diversity(metric, biom_data)
        alpha.append(c)
    # all_calcs = AlphaDiversityCalcs(calcs)
    # biom_data = all_calcs.getBiomData(biom)
    sample_metadata = dict(zip(biom_data.ids(), biom_data.metadata()))
    # result = all_calcs(data_path=biom_data)
    # return _format_alpha_results(result, sample_metadata=sample_metadata)
    return calcs


def alpha(table: biom.Table, metric: str):
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    counts = table.matrix_data.toarray().astype(float).T
    sample_ids = table.ids(axis='sample')

    result = skbio.diversity.alpha_diversity(metric=metric, counts=counts,
                                             ids=sample_ids)
    result.name = metric
    return result


print(alpha())


(array([[ 1.        ,  0.        , -0.        ],
       [ 1.        ,  0.        , -0.        ],
       [42.        ,  0.97297728,  5.24199479],
       [42.        ,  0.97297728,  5.24199479]]), array([u'test.fasta', u'GCF_000861845.1_ViralProj15432_genomic.fna',
       u'virus_setB_10K.fastq.fasta', 'virus_setB_10K.fastq.fasta (1)'],
      dtype=object), ['chao1', 'simpson', 'shannon'])

# {'labels': {u'label0': u'No Label'},
# 'data':
# [{'metadata': {'label': u'label0'}, 'values': {'shannon': -0.0, 'simpson': 0.0, 'chao1': 1.0}, 'sample_name': u'test.fasta'},
# {'metadata': {'label': u'label0'}, 'values': {'shannon': -0.0, 'simpson': 0.0, 'chao1': 1.0}, 'sample_name': u'GCF_000861845.1_ViralProj15432_genomic.fna'},
# {'metadata': {'label': u'label0'}, 'values': {'shannon': 5.241994787394689, 'simpson': 0.9729772752137023, 'chao1': 42.0}, 'sample_name': u'virus_setB_10K.fastq.fasta'},
# {'metadata': {'label': u'label0'}, 'values': {'shannon': 5.241994787394689, 'simpson': 0.9729772752137023, 'chao1': 42.0}, 'sample_name': 'virus_setB_10K.fastq.fasta (1)'}], 'methods': ['chao1', 'simpson', 'shannon']}

# def _format_alpha_results(result, sample_metadata=None):
#     if sample_metadata is None:
#         sample_metadata = {}
#
#     data, sample_names, calc_names = result
#     sample_list = []
#     labels_mapping = {}
#
#     for sample_index, sample_name in enumerate(sample_names):
#         single_sample_metadata = sample_metadata.get(sample_name, {})
#
#         sample_dict = dict(sample_name=sample_name, values=dict(),
#                            metadata=dict(label=single_sample_metadata.get('label', 'label0')))
#
#         labels_mapping[single_sample_metadata.get('label', 'label0')] = \
#             single_sample_metadata.get('label_name', 'No Label')
#
#         for calc_index, method_name in enumerate(calc_names):
#             value = data[sample_index][calc_index]
#
#             if str(value) == 'nan':
#                 # when sample has 0-hits, alpha diversity returns nan value
#                 value = 0
#
#             sample_dict['values'][method_name] = value
#
#         sample_list.append(sample_dict)
#
#     return dict(data=sample_list, methods=calc_names, labels=labels_mapping)