from biom import Table
from numpy import array

REGULAR_BIOM_TABLE = Table(data=array([
    [6.0, 0.0],
    [141.0, 67.0],
    [0.0, 6.0],
    [260.0, 601.0],
    [6128.0, 393.0],
    [35.0, 0.0],
    [0.0, 262.0],
    [0.0, 7.0],
    [19.0, 0.0]
]),
    observation_ids=['Dill cryptic virus 2', 'Enterobacteria phage T4', 'Hepatitis C virus',
                     'Human papillomavirus type 90', 'Lactobacillus phage Lv 1',
                     'Merkel cell polyomavirus', 'Mycobacterium phage Adler',
                     'Propionibacterium phge P105', 'Staphylococcus phage PH15'],
    sample_ids=['vag_intr_SRS014465.fasta', 'vag_intr_SRS015071.fasta'],
    sample_metadata=[{'name': 'vag_intr_SRS014465.fasta', 'file_id': '762b8657-bc6c-4c2f-8572-6be4df1adfc9',
                      'dataset_id': 'c1a84ab2-bca8-414b-9132-de4981426ba1', 'reads_total': 875954, 'label': 'label0',
                      'label_name': 'No Label'},
                     {'name': 'vag_intr_SRS015071.fasta', 'file_id': '8be2175c-6106-459f-859b-1e8f8bc6b0e8',
                      'dataset_id': '11d395b8-abfc-4571-aab6-e734b5d33885', 'reads_total': 507176, 'label': 'label0',
                      'label_name': 'No Label'}],
    table_id='05775205-8479-4346-9866-a45bdb449d70',
    type="OTU table")

REGULAR_BIOM_SAMPLE_META = {'vag_intr_SRS014465.fasta': {'id': 'vag_intr_SRS014465.fasta',
                                                         'metadata': {'name': 'vag_intr_SRS014465.fasta',
                                                                      'file_id': '762b8657-bc6c-4c2f-8572-6be4df1adfc9',
                                                                      'dataset_id': 'c1a84ab2-bca8-414b-9132-de4981426ba1',
                                                                      'reads_total': 875954, 'label': 'label0',
                                                                      'label_name': 'No Label'}},
                            'vag_intr_SRS015071.fasta': {'id': 'vag_intr_SRS015071.fasta',
                                                         'metadata': {'name': 'vag_intr_SRS015071.fasta',
                                                                      'file_id': '8be2175c-6106-459f-859b-1e8f8bc6b0e8',
                                                                      'dataset_id': '11d395b8-abfc-4571-aab6-e734b5d33885',
                                                                      'reads_total': 507176, 'label': 'label0',
                                                                      'label_name': 'No Label'}}}
