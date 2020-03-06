from os.path import isfile

import biom
from biom import load_table
from biom.table import Table
from numpy import asarray
from skbio.diversity import alpha_diversity

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


'''
BETA json output
{"nodes": {"GCF_000861845.1_ViralProj15432_genomic.fna": {"label": "label0"}, 
"GCF_000861845.1_ViralProj15432_genomic.fna (1)": {"label": "label0"}, 
"virus_setB_10K.fastq.fasta": {"label": "label0"}}, "labels": {"label0": "No Label"}, 
"data": [{"source": "GCF_000861845.1_ViralProj15432_genomic.fna (1)", "values": {"abund_jaccard": 1.0, "bray_curtis": 1.0}, "target": "virus_setB_10K.fastq.fasta"}, 
{"source": "GCF_000861845.1_ViralProj15432_genomic.fna", "values": {"abund_jaccard": 1.0, "bray_curtis": 1.0}, "target": "virus_setB_10K.fastq.fasta"}, 
{"source": "GCF_000861845.1_ViralProj15432_genomic.fna", "values": {"abund_jaccard": 0.0, "bray_curtis": 0.0}, "target": "GCF_000861845.1_ViralProj15432_genomic.fna (1)"}],
 "methods": ["abund_jaccard", "bray_curtis"]}
'''

'''
BETA LOGS 
 []
diversity.py:199 - {(u'virus_setB_10K.fastq.fasta', u'virus_setB_10K.fastq.fasta (1)'): {'abund_jaccard': 2.220446049250313e-16, 'bray_curtis': 0.0}, (u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta'): {'abund_jaccard': 1.0, 'bray_curtis': 1.0}, (u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta (1)'): {'abund_jaccard': 1.0, 'bray_curtis': 1.0}}
diversity.py:200 - {u'GCF_000861845.1_ViralProj15432_genomic.fna': {u'id': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'metadata': {u'name': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'reads_total': 1, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'4cc50293-3ac9-4565-98d2-0395f5ff14d1', u'dataset_id': u'e0d0e947-1d20-4ee0-918c-c763c429f8ad'}}, u'virus_setB_10K.fastq.fasta (1)': {u'id': u'virus_setB_10K.fastq.fasta (1)', u'metadata': {u'name': u'virus_setB_10K.fastq.fasta', u'reads_total': 390003, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'786b411b-c1a2-4142-87e3-e3f4286f72a0', u'dataset_id': u'e0f59313-ef3f-4b07-ae6e-a5377f610587'}}, u'virus_setB_10K.fastq.fasta': {u'id': u'virus_setB_10K.fastq.fasta', u'metadata': {u'name': u'virus_setB_10K.fastq.fasta', u'reads_total': 390003, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'0f0cbaa7-8f10-402c-9676-ec95048664f4', u'dataset_id': u'7d6e430c-85fa-4c6c-bd06-865c570feb8b'}}}
diversity.py:201 - {'abund_jaccard': [[0.0, 1.0, 1.0], [1.0, 0.0, 2.220446049250313e-16], [1.0, 2.220446049250313e-16, 0.0]], 'bray_curtis': [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]}
diversity.py:202 - [u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta', u'virus_setB_10K.fastq.fasta (1)']
'''

#BETA in quiime2
# def beta(table: biom.Table, metric: str,
#          pseudocount: int = 1, n_jobs: int = 1) -> skbio.DistanceMatrix:
#
#     if not (metric in non_phylogenetic_metrics()):
#         raise ValueError("Unknown metric: %s" % metric)
#
#     counts = table.matrix_data.toarray().T
#
#     def aitchison(x, y, **kwds):
#         return euclidean(clr(x), clr(y))
#
#     def canberra_adkins(x, y, **kwds):
#         if (x < 0).any() or (y < 0).any():
#             raise ValueError("Canberra-Adkins is only defined over positive "
#                              "values.")
#
#         nz = ((x > 0) | (y > 0))
#         x_ = x[nz]
#         y_ = y[nz]
#         nnz = nz.sum()
#
#         return (1. / nnz) * np.sum(np.abs(x_ - y_) / (x_ + y_))
#
#     if metric == 'aitchison':
#         counts += pseudocount
#         metric = aitchison
#     elif metric == 'canberra_adkins':
#         metric = canberra_adkins
#
#     if table.is_empty():
#         raise ValueError("The provided table object is empty")
#
#     sample_ids = table.ids(axis='sample')
#
#     return skbio.diversity.beta_diversity(
#         metric=metric,
#         counts=counts,
#         ids=sample_ids,
#         validate=True,
#         pairwise_func=sklearn.metrics.pairwise_distances,
#         n_jobs=n_jobs
#     )

analysis_alpha_diversity = alpha(input_biom)
print(analysis_alpha_diversity)
