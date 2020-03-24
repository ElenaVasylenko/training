import os
import logging
import json
import csv
import shutil
import unittest

MAPPING = {
    'phylogeny_analysis_tree.xml': 'phylogeny',
    'taxonomy_analysis_tree.xml': 'taxonomy',
}

LOGGER = logging.getLogger(__name__)
logging.basicConfig()


class

class BaseDiversityTestCase(unittest.TestCase):
    out_dir = 'tests/test_diversity/tmp_output'
    test_data_dir = 'tests/test_diversity/dummy_data'
    compare_data_dir = 'tests/test_diversity/to_compare'

    def setUp(self):
        os.makedirs(self.out_dir)

    def get_biom_file(self, biom_file):
        return os.path.join(self.test_data_dir, '{}.json'.format(biom_file))

    def get_tsv_file(self, tsv_name):
        return os.path.join(self.compare_data_dir, '{}.tsv'.format(tsv_name))

    def read_tsv_for_comparation(self, input_file):
        data = dict()
        with open(input_file, 'rb') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')
            col_names = tsv_reader.next()[1:]
            for row in tsv_reader:
                LOGGER.debug(row)
                row_name = row[0]
                values = row[1:]
                for value_index, value in enumerate(values):
                    data[(row_name, col_names[value_index])] = float(value) if value != 'nan' else 0
        return data

    def compare_alpha(self, json_out, original_output):
        tsv_data = self.read_tsv_for_comparation(original_output)
        with open(json_out) as file_:
            json_data = json.loads(file_.read())

            for value in json_data.get('data'):
                sample_name = value.get('sample_name')
                sample_values = value.get('values')

                for key, v in sample_values.items():
                    self.assertAlmostEquals(v, tsv_data.get((sample_name, key)), 4)

    def compare_beta(self, json_out, original_output_dir, file_name):
        tsv_data = {method: self.read_tsv_for_comparation(os.path.join(original_output_dir,
                                                                       '{}_{}.txt'.format(method, file_name)))
                    for method in ['abund_jaccard', 'bray_curtis']}

        with open(json_out) as file_:
            json_data = json.loads(file_.read())

            for value in json_data.get('data'):
                target = value.get('target')
                source = value.get('source')
                sample_values = value.get('values')

                for key, v in sample_values.items():
                    self.assertAlmostEquals(v, tsv_data[key][(target, source)], 4)

    def compare_pcoa(self, tsv_path, pcoa_path, method):
        tsv_data = self.read_tsv_for_comparation(tsv_path)

        with open(pcoa_path) as pcoa_file:
            data = json.loads(pcoa_file.read())
            samples_list = data.get('data')

        for sample in samples_list:
            sample_name = sample.get('sample_name')
            values_to_compare = sample.get('values')[method]

            for i, value in enumerate(values_to_compare):
                key = (str(sample_name), str(i))
                self.assertAlmostEquals(tsv_data[key], value, 4)

    def tearDown(self):
        shutil.rmtree(self.out_dir)
