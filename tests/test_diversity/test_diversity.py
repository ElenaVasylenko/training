import logging
import unittest

from diversity import alpha, BetaDiversityAnalysis
from .base import BaseDiversityTestCase

LOGGER = logging.getLogger(__name__)


class DiversityTestCase(BaseDiversityTestCase):

    def test_alpha_diversity(self):
        file_to_compare = 'alpha_regular_biom.json'
        biom_table = self.get_biom_table()
        alpha_diversity_result = alpha(biom_table)
        self.compare_alpha_sk_to_quime(alpha_diversity_result, file_to_compare)
        self.check_labels(alpha_diversity_result, file_to_compare)

    def test_beta_diversity(self):
        file_to_compare = 'beta_regular_biom.json'
        biom_table = self.get_biom_table()
        biom_meta = self.get_sample_meta()
        beta_diversity_analysis = BetaDiversityAnalysis(sample_meta_mapping=biom_meta)
        beta_diversity_object = beta_diversity_analysis.beta(table=biom_table)
        beta_diversity_json_result = beta_diversity_analysis.format_beta_diversity_to_json()
        self.compare_beta_sk_to_quime(beta_diversity_json_result, file_to_compare)
        self.check_labels(beta_diversity_json_result, file_to_compare)


if __name__ == '__main__':
    unittest.main()
