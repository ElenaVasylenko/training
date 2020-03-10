import json
from os.path import isfile

import biom
import skbio
import sklearn
from biom import load_table
from biom.table import Table
from numpy import asarray
from skbio import diversity
from skbio.diversity import alpha_diversity, beta_diversity

from diversity_data import input_biom

ALPHA_DIVERSITY_METHODS = ['chao1', 'simpson', 'shannon']
BETA_DIVERSITY_METHODS = ['jaccard', 'braycurtis']


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
    formatted_diversity_results = _format_alpha_results_to_json(aggregated_diversity_results, sample_metadata)

    return formatted_diversity_results


def aggregate_results(table, sample_ids):
    """

    :param table:
    :param sample_ids:
    :return:
    """
    diversity_aggregated = (asarray(table).transpose(), sample_ids, ALPHA_DIVERSITY_METHODS)
    return diversity_aggregated


def _format_alpha_results_to_json(result, sample_metadata=None):
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
data                                   []
data_mapping        diversity.py:199 - {(u'virus_setB_10K.fastq.fasta', u'virus_setB_10K.fastq.fasta (1)'): 
{'abund_jaccard': 2.220446049250313e-16, 'bray_curtis': 0.0}, (u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta'): {'abund_jaccard': 1.0, 'bray_curtis': 1.0}, (u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta (1)'): {'abund_jaccard': 1.0, 'bray_curtis': 1.0}}
sample_meta_mapping diversity.py:200 - {u'GCF_000861845.1_ViralProj15432_genomic.fna': {u'id': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'metadata': {u'name': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'reads_total': 1, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'4cc50293-3ac9-4565-98d2-0395f5ff14d1', u'dataset_id': u'e0d0e947-1d20-4ee0-918c-c763c429f8ad'}}, u'virus_setB_10K.fastq.fasta (1)': {u'id': u'virus_setB_10K.fastq.fasta (1)', u'metadata': {u'name': u'virus_setB_10K.fastq.fasta', u'reads_total': 390003, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'786b411b-c1a2-4142-87e3-e3f4286f72a0', u'dataset_id': u'e0f59313-ef3f-4b07-ae6e-a5377f610587'}}, u'virus_setB_10K.fastq.fasta': {u'id': u'virus_setB_10K.fastq.fasta', u'metadata': {u'name': u'virus_setB_10K.fastq.fasta', u'reads_total': 390003, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'0f0cbaa7-8f10-402c-9676-ec95048664f4', u'dataset_id': u'7d6e430c-85fa-4c6c-bd06-865c570feb8b'}}}
raw_data            diversity.py:201 - {'abund_jaccard': [[0.0, 1.0, 1.0], [1.0, 0.0, 2.220446049250313e-16], [1.0, 2.220446049250313e-16, 0.0]], 
                                         'bray_curtis': [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]}
raw_ids             diversity.py:202 - [u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta', u'virus_setB_10K.fastq.fasta (1)']
'''

#
# def beta(table: biom.Table):
#     if table.is_empty():
#         raise ValueError("The provided table object is empty")
#
#     beta_diversities = []
#     counts = table.matrix_data.toarray().astype(float).T
#     print(counts)
#     sample_ids = table.ids(axis='sample')
#     print(sample_ids)
#     sample_metadata = dict(zip(table.ids(), table.metadata()))
#
#     for metric in BETA_DIVERSITY_METHODS:
#         result = beta_diversity(metric=metric, counts=counts, ids=sample_ids, pairwise_func=sklearn.metrics.pairwise_distances)
#         result.name = metric
#         print(result.data)
#         beta_diversities.append(result)
#     print(beta_diversities)
#     return beta_diversities
#
#
# beta_div = beta(input_biom)
# print(beta_div)
# print(diversity.get_beta_diversity_metrics())
# for res in beta_div:
#     print(res.data)
#     print(res.ids)
#     print(res.name)

# analysis_alpha_diversity = alpha(input_biom)
# print(analysis_alpha_diversity)

# BETA in quiime2
#
#
# def write_data(data):
#     data, ids = data
#     method = self.methods[self.method_index]
#     for index_x, id_from in enumerate(ids):
#         for index_y, id_to in enumerate(ids):
#             if index_y <= index_x:
#                 continue
#
#             key = (id_from, id_to)
#             if key not in self.data_mapping:
#                 self.data_mapping[key] = dict()
#             self.data_mapping[key][method] = data[index_x][index_y]
#
#     self.raw_data[method] = data
#     self.raw_ids = ids
#     self.method_index += 1
#
# def flush(self):
#     data = []
#     for key, methods in self.data_mapping.items():
#         source, target = key
#         link = dict(source=source, target=target, values=dict())
#         for method, value in methods.items():
#             link['values'][method] = value
#
#         data.append(link)
#
#     with open(self.target_file, 'w') as file_:
#         labels_mapping = {sample['metadata']['label']: sample['metadata']['label_name']
#                           for sample in self.sample_meta_mapping.values()}
#         nodes_info = {sample_id: dict(label=sample['metadata']['label'])
#                       for sample_id, sample in self.sample_meta_mapping.items()}
#
#         file_.write_data(json.dumps(dict(data=data, methods=self.methods, labels=labels_mapping, nodes=nodes_info)))
#
#     return self.data_mapping

methods = ['jaccard', 'braycurtis']

data_mapping = {}


def beta(table, metric: str, pseudocount: int = 1, n_jobs: int = 1):
    counts = table.matrix_data.toarray().T

    if table.is_empty():
        raise ValueError("The provided table object is empty")

    sample_ids = table.ids(axis='sample')
    # for metric in BETA_DIVERSITY_METHODS:
    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        validate=True,
        pairwise_func=sklearn.metrics.pairwise_distances,
        n_jobs=n_jobs
    )


def format_beta_diversity_distance_matrix(ids, data):
    return list(data), list(ids)


def get_beta_diversity_for_pair_files_intersection(data, methods):
    data, ids = data.data, data.ids

    for method in methods:
        for index_x, id_from in enumerate(ids):
            for index_y, id_to in enumerate(ids):

                if index_y <= index_x:
                    continue

                key = (id_from, id_to)
                if key not in data_mapping:
                    data_mapping[key] = dict()
                data_mapping[key][method] = data[index_x][index_y]

    return data_mapping


def get_sample_meta_mapping(biom_file):
    with open(biom_file) as file_:
        sample_meta_mapping = {sample['id']: sample for sample in json.loads(file_.read()).get('columns')}
    return sample_meta_mapping


def format_beta_diversity_to_json(data_mapping_items, sample_meta_mapping, target_file):
    data = []
    for key, methods in data_mapping.items():
        source, target = key
        link = dict(source=source, target=target, values=dict())
        for method, value in methods.items():
            link['values'][method] = value

        data.append(link)

    with open(target_file, 'w') as file_:
        labels_mapping = {sample['metadata']['label']: sample['metadata']['label_name']
                          for sample in sample_meta_mapping.values()}
        nodes_info = {sample_id: dict(label=sample['metadata']['label'])
                      for sample_id, sample in sample_meta_mapping.items()}


        file_.write(json.dumps(dict(data=data, methods=BETA_DIVERSITY_METHODS, labels=labels_mapping, nodes=nodes_info)))

    return data_mapping

# b = beta(input_biom, metric='braycurtis')
b2 = beta(input_biom, metric='jaccard')
# print(b.data)
print(b2.data)
# print(format_beta_diversity_distance_matrix(b.ids, b.data))
print(format_beta_diversity_distance_matrix(b2.ids, b2.data))

beta_data = get_beta_diversity_for_pair_files_intersection(b2, methods)
print(beta_data)


#[[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0]],
# [u'GCF_000861845.1_ViralProj15432_genomic.fna', u'GCF_000861845.1_ViralProj15432_genomic.fna (1)', u'virus_setB_10K.fastq.fasta'])
# {"nodes": {"GCF_000861845.1_ViralProj15432_genomic.fna": {"label": "label0"},
#            "GCF_000861845.1_ViralProj15432_genomic.fna (1)": {"label": "label0"},
#            "virus_setB_10K.fastq.fasta": {"label": "label0"}}, "labels": {"label0": "No Label"},
#  "data": [{"source": "GCF_000861845.1_ViralProj15432_genomic.fna (1)", "values": {"abund_jaccard": 1.0, "bray_curtis": 1.0}, "target": "virus_setB_10K.fastq.fasta"},
#           {"source": "GCF_000861845.1_ViralProj15432_genomic.fna", "values": {"abund_jaccard": 1.0, "bray_curtis": 1.0}, "target": "virus_setB_10K.fastq.fasta"},
#           {"source": "GCF_000861845.1_ViralProj15432_genomic.fna", "values": {"abund_jaccard": 0.0, "bray_curtis": 0.0}, "target": "GCF_000861845.1_ViralProj15432_genomic.fna (1)"}
#          ], "methods": ["abund_jaccard", "bray_curtis"]}