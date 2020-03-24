import os
import json
import logging
from diversity import alpha, BetaDiversityAnalysis
from .base import BaseDiversityTestCase

LOGGER = logging.getLogger(__name__)
from comparative_helper import ComparativeAnalysis


#
def alpha_diversity(self):
    matrix_data = self._build_data_matrix_for_field(FieldsForComparative.TOTAL_HIT_FREQUENCY,
                                                    dtype=np.float64,
                                                    lg=False)
    biom_table = self.get_biom_table(self._init_table(matrix_data))
    out_path = os.path.join(self._workdir, 'alpha.json')
    self.logger.debug("Path to store alpha %s, exclude path %s", out_path, self.exclude_prefix)

    alpha_diversity_result = alpha(biom_table)
    with open(out_path, 'w') as file_:
        file_.write(json.dumps(alpha_diversity_result))

    return os.path.relpath(out_path, self.exclude_prefix)


def beta_diversity(self):
    biom_file = self.biom_v1_path
    out_path = os.path.join(self._workdir, 'beta.json')
    matrix_data = self._build_data_matrix_for_field(FieldsForComparative.TOTAL_HIT_FREQUENCY,
                                                    dtype=np.float64,
                                                    lg=False)
    biom_table = self.get_biom_table(self._init_table(matrix_data))

    with open(biom_file) as file_:
        sample_meta_mapping = {sample['id']: sample for sample in json.loads(file_.read()).get('columns')}

    beta_diversity_analysis = BetaDiversityAnalysis(sample_meta_mapping=sample_meta_mapping)
    self.beta_diversity_object = beta_diversity_analysis.beta(table=biom_table)
    beta_diversity_json_result = beta_diversity_analysis.format_beta_diversity_to_json()
    with open(out_path, 'w') as file_:
        file_.write(json.dumps(beta_diversity_json_result))

    return os.path.relpath(out_path, self.exclude_prefix)


class DiversityTestCase(BaseDiversityTestCase):


    def test_0_hit_biom_alpha(self):
        """
        Test generating alpha diversity for biom with 0-hit analysis.

        Input - 0-hit biom
        """
        input_biom = self.get_biom_file('0-hit-biom')
        file_to_compare = '0-hit-biom'
        alpha_diversity_result = alpha(input_biom)
        self.compare_alpha(alpha_diversity_result, file_to_compare)

    def test_0_hit_biom_beta(self):
        """
        Test monkey patching for beta works as expected with 0-hit-biom.

        Input - Usual biom
        Expected - tsv and json data should be logicaly identical
        """
        input_biom = self.get_biom_file('0-hit-biom')
        output_json = os.path.join(self.out_dir, 'beta.json')

        beta_diversity_analysis = BetaDiversityAnalysis(sample_meta_mapping=sample_meta_mapping)
        self.beta_diversity_object = beta_diversity_analysis.beta(table=biom_table)
        beta_diversity_json_result = beta_diversity_analysis.format_beta_diversity_to_json()

        self.compare_beta(output_json, self.out_dir, '0-hit-biom')

    def test_0_hit_two_samples_beta(self):
        """
        Test monkey patching for beta works as expected with 0-hit-biom, where 2 samples have 0 hits.

        Input - Usual biom
        Expected - tsv and json data should be logicaly identical
        """
        input_biom = self.get_biom_file('0-hit-3-samples-biom')
        output_json = os.path.join(self.out_dir, 'beta.json')

        beta_diversity(input_biom, output_json)

        self.compare_beta(output_json, self.out_dir, '0-hit-3-samples-biom')

    def test_beta_labels(self):
        """
        Test if labels for samples in beta diversity matrix are assigned correctly.

        Input - Usual biom
        Expected - label id in beta matrix should be the same as in biom
        """
        input_biom = self.get_biom_file('common_biom')
        output_json = os.path.join(self.out_dir, 'beta.json')

        beta_diversity(input_biom, output_json)
        # single_file_beta(input_biom, ['abund_jaccard', 'bray_curtis'], None, self.out_dir, None)

        with open(input_biom) as file_:
            input_biom_data = json.loads(file_.read())

        with open(output_json) as file_:
            beta_raw = file_.read()
            LOGGER.debug(beta_raw)
            output_beta_matrix = json.loads(beta_raw)
            output_beta_nodes = output_beta_matrix.get('nodes')
            output_beta_labels = output_beta_matrix.get('labels')

        sample_id_to_sample = {sample['id']: sample.get('metadata') for sample in input_biom_data.get('columns')}

        for sample_id, sample_data in output_beta_nodes.items():
            label_id = sample_data.get('label')
            self.assertEquals(sample_id_to_sample.get(sample_id).get('label'), label_id)
            self.assertEquals(sample_id_to_sample.get(sample_id).get('label_name'), output_beta_labels.get(label_id))

    def test_alpha_labels(self):
        """
        Test if labels for samples in alpha diversity are assigned correctly.

        Input - Usual biom
        Expected - label id in alpha diversity should be the same as in biom
        """
        input_biom = self.get_biom_file('common_biom')
        output_json = os.path.join(self.out_dir, 'alpha.json')

        alpha_diversity(input_biom, output_json)

        with open(input_biom) as file_:
            input_biom_data = json.loads(file_.read())

        with open(output_json) as file_:
            alpha_raw = file_.read()
            LOGGER.debug(alpha_raw)
            output_alpha = json.loads(alpha_raw)
            output_alpha_nodes = output_alpha.get('data')
            output_alpha_labels = output_alpha.get('labels')

        sample_id_to_sample = {sample['id']: sample.get('metadata') for sample in input_biom_data.get('columns')}

        for sample_data in output_alpha_nodes:
            label_id = sample_data['metadata']['label']
            sample_id = sample_data['sample_name']
            self.assertEquals(sample_id_to_sample.get(sample_id).get('label'), label_id)
            self.assertEquals(sample_id_to_sample.get(sample_id).get('label_name'), output_alpha_labels.get(label_id))
