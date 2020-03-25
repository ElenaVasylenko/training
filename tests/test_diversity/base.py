import json
import os
import unittest
from enum import Enum

from .table_data.tables import REGULAR_BIOM_SAMPLE_META, REGULAR_BIOM_TABLE


class DiversityResult(Enum):
    DATA = 'data'
    SAMPLE_NAME = 'sample_name'
    VALUES = 'values'
    TARGET = 'target'
    SOURCE = 'source'
    LABELS = 'labels'
    NODES = 'nodes'


class BaseDiversityTestCase(unittest.TestCase):
    compare_results_dir = '/home/ovasylenko/PycharmProjects/training/tests/test_diversity/quiime1_result_data'
    # compare_results_dir = 'tests/test_diversity/quiime1_result_data'
    ALPHA_PRECISION_PARAMETER = 2
    BETA_PRECISION_PARAMETER = 1

    def get_biom_table(self):
        return REGULAR_BIOM_TABLE

    def get_sample_meta(self):
        return REGULAR_BIOM_SAMPLE_META

    def compare_alpha_sk_to_quime(self, skbio_output, quime_output):
        quime_output = os.path.join(self.compare_results_dir, quime_output)

        with open(quime_output) as file_:
            json_data = json.loads(file_.read())
            for value in json_data.get(DiversityResult.DATA.value):
                for sk_value in skbio_output.get(DiversityResult.DATA.value):
                    if value.get(DiversityResult.SAMPLE_NAME.value) == sk_value.get(DiversityResult.SAMPLE_NAME.value):
                        sample_values = value.get(DiversityResult.VALUES.value)
                        skbio_sample_values = sk_value.get(DiversityResult.VALUES.value)
                        break
                for key, v in sample_values.items():
                    self.assertAlmostEqual(v, skbio_sample_values[key], self.ALPHA_PRECISION_PARAMETER)

    def compare_beta_sk_to_quime(self, skbio_output, quime_output):
        quime_output = os.path.join(self.compare_results_dir, quime_output)

        with open(quime_output) as file_:
            json_data = json.loads(file_.read())
            for value in json_data.get(DiversityResult.DATA.value):
                for sk_value in skbio_output.get(DiversityResult.DATA.value):
                    if value.get(DiversityResult.TARGET.value) == sk_value.get(DiversityResult.TARGET.value) \
                            and value.get(DiversityResult.SOURCE.value) == sk_value.get(DiversityResult.SOURCE.value):
                        skbio_sample_values = sk_value.get(DiversityResult.VALUES.value)
                        sample_values = value.get(DiversityResult.VALUES.value)
                        break
                for key, v in sample_values.items():
                    self.assertAlmostEqual(skbio_sample_values[key], v, delta=self.BETA_PRECISION_PARAMETER)

    def check_labels(self, skbio_output, quime_output):
        """
        Test if labels for samples in beta diversity matrix are assigned correctly.

        Input - Usual biom
        Expected - label id in beta matrix should be the same as in biom
        """
        quime_output = os.path.join(self.compare_results_dir, quime_output)

        with open(quime_output) as file_:
            beta_raw = file_.read()
            output_beta_matrix = json.loads(beta_raw)
            output_beta_nodes = output_beta_matrix.get(DiversityResult.NODES.value)
            output_beta_labels = output_beta_matrix.get(DiversityResult.LABELS.value)
        sk_labels = skbio_output.get(DiversityResult.LABELS.value)
        sk_nodes = skbio_output.get(DiversityResult.NODES.value)
        assert output_beta_labels == sk_labels
        assert output_beta_nodes == sk_nodes
