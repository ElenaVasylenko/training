import json
import biom
import sklearn
from skbio.diversity import alpha_diversity, beta_diversity
from diversity_data import input_biom

#
# ALPHA_DIVERSITY_METHODS = ['chao1', 'simpson', 'shannon']
BETA_DIVERSITY_METHODS = ['jaccard', 'braycurtis']


class BetaDiversityAnalysis(object):

    def __init__(self, target_file, methods, biom_file):
        self.methods = methods
        self.method_index = 0
        self.data = []
        self.data_mapping = dict()
        self.target_file = target_file
        self.raw_ids = None
        self.raw_data = {}

        self.sample_meta_mapping = {u'GCF_000861845.1_ViralProj15432_genomic.fna': {u'id': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'metadata': {u'name': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'reads_total': 1, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'4cc50293-3ac9-4565-98d2-0395f5ff14d1', u'dataset_id': u'e0d0e947-1d20-4ee0-918c-c763c429f8ad'}}, u'virus_setB_10K.fastq.fasta': {u'id': u'virus_setB_10K.fastq.fasta', u'metadata': {u'name': u'virus_setB_10K.fastq.fasta', u'reads_total': 390003, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'1af0ba9b-b94a-485e-9232-37dffc48754f', u'dataset_id': u'db1518c1-8779-423b-9f48-92d1c2bb77a3'}}, u'GCF_000861845.1_ViralProj15432_genomic.fna (1)': {u'id': u'GCF_000861845.1_ViralProj15432_genomic.fna (1)', u'metadata': {u'name': u'GCF_000861845.1_ViralProj15432_genomic.fna', u'reads_total': 1, u'label': u'label0', u'label_name': u'No Label', u'file_id': u'29c6d2e6-617c-494c-b868-6b6ec5518d81', u'dataset_id': u'9c27b10c-687e-43d2-9631-21a416489541'}}}


    def compute_beta_diversity(self):
        for metric in self.methods:
            beta_diversity_result = self.beta(input_biom, metric=metric)
            self.get_beta_diversity_for_pair_files_intersection(beta_diversity_result)
        self.format_beta_diversity_to_json()

    def beta(self, table: biom.Table, metric: str, pseudocount: int = 1, n_jobs: int = 1):
        counts = table.matrix_data.toarray().T

        if table.is_empty():
            raise ValueError("The provided table object is empty")

        sample_ids = table.ids(axis='sample')
        beta_dv = beta_diversity(
            metric=metric,
            counts=counts,
            ids=sample_ids,
            validate=True,
            pairwise_func=sklearn.metrics.pairwise_distances,
            n_jobs=n_jobs
        )
        return beta_dv

    def get_beta_diversity_for_pair_files_intersection(self, data):

        data, ids = data.data, data.ids
        method = self.methods[self.method_index]
        for index_x, id_from in enumerate(ids):
            for index_y, id_to in enumerate(ids):
                if index_y <= index_x:
                    continue

                key = (id_from, id_to)
                if key not in self.data_mapping:
                    self.data_mapping[key] = dict()
                self.data_mapping[key][method] = data[index_x][index_y]

        self.raw_data[method] = data
        self.raw_ids = ids
        self.method_index += 1

    def get_all_metrics(self):
        return self.raw_ids, self.raw_data

    def format_beta_diversity_to_json(self):
        data = []
        for key, methods in self.data_mapping.items():
            source, target = key
            link = dict(source=source, target=target, values=dict())
            for method, value in methods.items():
                link['values'][method] = value

            data.append(link)

        with open(self.target_file, 'w') as file_:
            labels_mapping = {sample['metadata']['label']: sample['metadata']['label_name']
                              for sample in self.sample_meta_mapping.values()}
            nodes_info = {sample_id: dict(label=sample['metadata']['label'])
                          for sample_id, sample in self.sample_meta_mapping.items()}
            print('lol')
            print(data, '\n')
            print(BETA_DIVERSITY_METHODS, '\n')
            print(labels_mapping, '\n')
            print(nodes_info, '\n')
            file_.write(json.dumps(dict(data=data, methods=self.methods, labels=labels_mapping, nodes=nodes_info)))
        print(json.dumps(dict(data=data, methods=self.methods, labels=labels_mapping, nodes=nodes_info)))
        return self.data_mapping


def monkey_patch_format_matrix(ids, data):
    return data.tolist(), ids.tolist()

beta_diversity_analysis = BetaDiversityAnalysis('result2.json', BETA_DIVERSITY_METHODS, input_biom)
beta_diversity_analysis.compute_beta_diversity()