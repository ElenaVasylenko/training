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
observation_ids = ['Abelson murine leukemia virus GCF 000848265.1',
                   'Avian orthoreovirus AVS B GCF 000891595.1',
                   'Camelpox virus GCF 000839105.1',
                   'Canine morbillivirus GCF 000854065.1',
                   'Colorado tick fever virus GCF 000853265.1',
                   'Cowpox virus Brighton Red GCF 000839185.1',
                   'Crimean Congo hemorrhagic fever orthonairovirus GCF 000854165.1',
                   'Cynomolgus cytomegalovirus GCF 001963235.1',
                   'Hepacivirus C GCF 000861845.1',
                   'Human T lymphotropic virus 2 GCF 000847505.1',
                   'Human T lymphotropic virus 4 GCF 000882595.1',
                   'Human bocavirus 2c PK GCF 000882675.1',
                   'Human mastadenovirus A 6017 Branch',
                   'Human metapneumovirus GCF 000865625.1',
                   'Human orthopneumovirus B1 GCF 000855545.1',
                   'Human papillomavirus 4 GCF 000864845.1',
                   'Human parvovirus B19 GCF 000839645.1',
                   'Human respirovirus 1 GCF 000848705.1',
                   'Influenza A virus A Puerto Rico 8 1934H1N1 GCF 000865725.1',
                   'JC polyomavirus Mad1 GCF 000863805.1',
                   'Lassa mammarenavirus Josiah GCF 000851705.1',
                   'Long Island tick rhabdovirus GCF 000926995.1',
                   'Lymphocytic choriomeningitis mammarenavirus GCF 000851025.1',
                   'Marburg marburgvirus GCF 000857325.2',
                   'Measles morbillivirus Ichinose B95a GCF 000854845.1',
                   'Menghai flavivirus GCF 002029595.1',
                   'Molluscum contagiosum virus subtype 1 GCF 000843325.1',
                   'Monkeypox virus Zaire 96 I 16 GCF 000857045.1',
                   'Mumps rubulavirus Miyahara GCF 000856685.1',
                   'Orbivirus SX 2017a GCF 002005745.1',
                   'Orf virus OV SA00 GCF 000844845.1',
                   'Pigeon picornavirus B GCF 000892795.1',
                   'Rabies lyssavirus GCF 000859625.1',
                   'Respiratory syncytial virus GCF 000856445.1',
                   'Rift Valley fever virus GCF 000847345.1',
                   'Rubella virus F Therien GCF 000863025.1',
                   'Saint Louis encephalitis virus GCF 000866785.1',
                   'Tick borne encephalitis virus GCF 000863125.1',
                   'Variola virus India 1967 ssp major GCF 000859885.1',
                   'Venezuelan equine encephalitis virus GCF 000862105.1',
                   'Wenling toga like virus GCF 001959035.1',
                   'Western equine encephalitis virus GCF 000850885.1',
                   'Whitewater Arroyo mammarenavirus GCF 000872865.1']
data = array([[0.,     0.,        5.29367628],
 [0.,         0.,         5.52763892],
 [0.,        0.,         1.43136376],
 [0.,         0.,         5.55793407],
 [0. ,        0.  ,       5.53187088],
 [0.  ,       0.  ,       1.69897   ],
 [0.   ,      0.     ,    5.56313327],
 [0.    ,     0.     ,    5.11045839],
 [3.53096768 , 3.53096768 , 0.],
 [0. ,        0. ,        5.48072107],
 [0.  ,       0.   ,      5.45509927],
 [0.     ,    0.  ,       5.50020628],
 [0.  ,       0. ,        5.50599603],
 [0.     ,    0.     ,    5.59733532],
 [0. ,        0.  ,       1.11394335],
 [0.  ,       0.   ,      5.58181886],
 [0.  ,       0. ,        5.51705331],
 [0.  ,       0. ,        5.42805538],
 [0.  ,       0.  ,       5.49400981],
 [0.   ,      0.   ,      5.56329831],
 [0.   ,      0.  ,       5.41639084],
 [0. ,        0.   ,      5.5121452 ],
 [0.  ,       0.  ,       5.56240945],
 [0.   ,      0.  ,       5.58252697],
 [0.   ,      0.  ,       5.52382058],
 [0.  ,       0.      ,   5.50203362],
 [0.  ,       0.    ,     5.39554654],
 [0.  ,       0.  ,       5.1721766 ],
 [0.  ,       0.     ,    5.55956312],
 [0.     ,    0.  ,       5.42083006],
 [0.   ,      0.  ,       5.29768831],
 [0.   ,      0.  ,       5.53319355],
 [0.   ,      0.    ,     5.54272586],
 [0.  ,       0.     ,    5.61124096],
 [0.    ,     0.  ,       5.54092606],
 [0.    ,     0.   ,      5.31292593],
 [0.  ,       0.    ,     5.50991456],
 [0.    ,     0.  ,       5.44642842],
 [0.  ,       0.    ,     5.12015589],
 [0.    ,     0.  ,       5.50973867],
 [0.  ,       0.       ,  5.4030879 ],
 [0.  ,       0.  ,       5.51752842],
 [0.   ,      0.   ,      5.57927586]])

input_biom = Table(data=array([[0.,     0.,        5.29367628],
 [0.,         0.,         5.52763892],
 [0.,        0.,         1.43136376],
 [0.,         0.,         5.55793407],
 [0. ,        0.  ,       5.53187088],
 [0.  ,       0.  ,       1.69897   ],
 [0.   ,      0.     ,    5.56313327],
 [0.    ,     0.     ,    5.11045839],
 [3.53096768 , 3.53096768 , 0.],
 [0. ,        0. ,        5.48072107],
 [0.  ,       0.   ,      5.45509927],
 [0.     ,    0.  ,       5.50020628],
 [0.  ,       0. ,        5.50599603],
 [0.     ,    0.     ,    5.59733532],
 [0. ,        0.  ,       1.11394335],
 [0.  ,       0.   ,      5.58181886],
 [0.  ,       0. ,        5.51705331],
 [0.  ,       0. ,        5.42805538],
 [0.  ,       0.  ,       5.49400981],
 [0.   ,      0.   ,      5.56329831],
 [0.   ,      0.  ,       5.41639084],
 [0. ,        0.   ,      5.5121452 ],
 [0.  ,       0.  ,       5.56240945],
 [0.   ,      0.  ,       5.58252697],
 [0.   ,      0.  ,       5.52382058],
 [0.  ,       0.      ,   5.50203362],
 [0.  ,       0.    ,     5.39554654],
 [0.  ,       0.  ,       5.1721766 ],
 [0.  ,       0.     ,    5.55956312],
 [0.     ,    0.  ,       5.42083006],
 [0.   ,      0.  ,       5.29768831],
 [0.   ,      0.  ,       5.53319355],
 [0.   ,      0.    ,     5.54272586],
 [0.  ,       0.     ,    5.61124096],
 [0.    ,     0.  ,       5.54092606],
 [0.    ,     0.   ,      5.31292593],
 [0.  ,       0.    ,     5.50991456],
 [0.    ,     0.  ,       5.44642842],
 [0.  ,       0.    ,     5.12015589],
 [0.    ,     0.  ,       5.50973867],
 [0.  ,       0.       ,  5.4030879 ],
 [0.  ,       0.  ,       5.51752842],
 [0.   ,      0.   ,      5.57927586]]),
observation_ids = ['Abelson murine leukemia virus GCF 000848265.1',
                   'Avian orthoreovirus AVS B GCF 000891595.1',
                   'Camelpox virus GCF 000839105.1',
                   'Canine morbillivirus GCF 000854065.1',
                   'Colorado tick fever virus GCF 000853265.1',
                   'Cowpox virus Brighton Red GCF 000839185.1',
                   'Crimean Congo hemorrhagic fever orthonairovirus GCF 000854165.1',
                   'Cynomolgus cytomegalovirus GCF 001963235.1',
                   'Hepacivirus C GCF 000861845.1',
                   'Human T lymphotropic virus 2 GCF 000847505.1',
                   'Human T lymphotropic virus 4 GCF 000882595.1',
                   'Human bocavirus 2c PK GCF 000882675.1',
                   'Human mastadenovirus A 6017 Branch',
                   'Human metapneumovirus GCF 000865625.1',
                   'Human orthopneumovirus B1 GCF 000855545.1',
                   'Human papillomavirus 4 GCF 000864845.1',
                   'Human parvovirus B19 GCF 000839645.1',
                   'Human respirovirus 1 GCF 000848705.1',
                   'Influenza A virus A Puerto Rico 8 1934H1N1 GCF 000865725.1',
                   'JC polyomavirus Mad1 GCF 000863805.1',
                   'Lassa mammarenavirus Josiah GCF 000851705.1',
                   'Long Island tick rhabdovirus GCF 000926995.1',
                   'Lymphocytic choriomeningitis mammarenavirus GCF 000851025.1',
                   'Marburg marburgvirus GCF 000857325.2',
                   'Measles morbillivirus Ichinose B95a GCF 000854845.1',
                   'Menghai flavivirus GCF 002029595.1',
                   'Molluscum contagiosum virus subtype 1 GCF 000843325.1',
                   'Monkeypox virus Zaire 96 I 16 GCF 000857045.1',
                   'Mumps rubulavirus Miyahara GCF 000856685.1',
                   'Orbivirus SX 2017a GCF 002005745.1',
                   'Orf virus OV SA00 GCF 000844845.1',
                   'Pigeon picornavirus B GCF 000892795.1',
                   'Rabies lyssavirus GCF 000859625.1',
                   'Respiratory syncytial virus GCF 000856445.1',
                   'Rift Valley fever virus GCF 000847345.1',
                   'Rubella virus F Therien GCF 000863025.1',
                   'Saint Louis encephalitis virus GCF 000866785.1',
                   'Tick borne encephalitis virus GCF 000863125.1',
                   'Variola virus India 1967 ssp major GCF 000859885.1',
                   'Venezuelan equine encephalitis virus GCF 000862105.1',
                   'Wenling toga like virus GCF 001959035.1',
                   'Western equine encephalitis virus GCF 000850885.1',
                   'Whitewater Arroyo mammarenavirus GCF 000872865.1'],
sample_ids = ['GCF_000861845.1_ViralProj15432_genomic.fna', 'GCF_000861845.1_ViralProj15432_genomic.fna (1)', 'virus_setB_10K.fastq.fasta'],
observation_metadata=[{'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Ortervirales', 'Retroviridae', 'Gammaretrovirus', 'Abelson murine leukemia virus', 'Abelson murine leukemia virus_u_t'], 'tax_id': 11788, 'title': 'Abelson murine leukemia virus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Reoviridae', 'Orthoreovirus', 'Avian orthoreovirus', 'Avian orthoreovirus'], 'tax_id': 38170, 'title': 'Avian orthoreovirus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Poxviridae', 'Orthopoxvirus', 'Camelpox virus', 'Camelpox virus'], 'tax_id': 28873, 'title': 'Camelpox virus'},
                      {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Paramyxoviridae', 'Morbillivirus', 'Canine morbillivirus', 'Canine morbillivirus'], 'tax_id': 11232, 'title': 'Canine morbillivirus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Reoviridae', 'Coltivirus', 'Colorado tick fever virus', 'Colorado tick fever virus'], 'tax_id': 46839, 'title': 'Colorado tick fever virus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Poxviridae', 'Orthopoxvirus', 'Cowpox virus', 'Cowpox virus'], 'tax_id': 10243, 'title': 'Cowpox virus'},
                      {'taxonomy': ['Viruses', 'Negarnaviricota', 'Ellioviricetes', 'Bunyavirales', 'Nairoviridae', 'Orthonairovirus', 'Crimean-Congo hemorrhagic fever orthonairovirus', 'Crimean-Congo hemorrhagic fever orthonairovirus'], 'tax_id': 1980519, 'title': 'Crimean-Congo hemorrhagic fever orthonairovirus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Herpesvirales', 'Herpesviridae', 'Cytomegalovirus', 'Cynomolgus cytomegalovirus', 'Cynomolgus cytomegalovirus'], 'tax_id': 1919083, 'title': 'Cynomolgus cytomegalovirus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Flaviviridae', 'Hepacivirus', 'Hepacivirus C', 'Hepatitis C virus genotype 1'], 'tax_id': 41856, 'title': 'Hepatitis C virus genotype 1'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Ortervirales', 'Retroviridae', 'Deltaretrovirus', 'Primate T-lymphotropic virus 2', 'Human T-lymphotropic virus 2'], 'tax_id': 11909, 'title': 'Human T-lymphotropic virus 2'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Ortervirales', 'Retroviridae', 'Deltaretrovirus', 'Human T-lymphotropic virus 4', 'Human T-lymphotropic virus 4'], 'tax_id': 318279, 'title': 'Human T-lymphotropic virus 4'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Parvoviridae', 'Bocaparvovirus', 'Primate bocaparvovirus 2', 'Human bocavirus 2c PK'], 'tax_id': 1511882, 'title': 'Human bocavirus 2c PK'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Adenoviridae', 'Mastadenovirus', 'Human mastadenovirus A', 'Human mastadenovirus A'], 'tax_id': 129875, 'title': 'Human mastadenovirus A'},
                      {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Pneumoviridae', 'Metapneumovirus', 'Human metapneumovirus', 'Human metapneumovirus'], 'tax_id': 162145, 'title': 'Human metapneumovirus'},
                      {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Pneumoviridae', 'Orthopneumovirus', 'Human orthopneumovirus', 'Human orthopneumovirus'], 'tax_id': 11250, 'title': 'Human orthopneumovirus'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Papillomaviridae', 'Gammapapillomavirus', 'Gammapapillomavirus 1', 'Human papillomavirus 4'], 'tax_id': 10617, 'title': 'Human papillomavirus 4'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Parvoviridae', 'Erythroparvovirus', 'Primate erythroparvovirus 1', 'Human parvovirus B19'], 'tax_id': 10798, 'title': 'Human parvovirus B19'},
                      {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Paramyxoviridae', 'Respirovirus', 'Human respirovirus 1', 'Human respirovirus 1'], 'tax_id': 12730, 'title': 'Human respirovirus 1'},
                      {'taxonomy': ['Viruses', 'Negarnaviricota', 'Insthoviricetes', 'Articulavirales', 'Orthomyxoviridae', 'Alphainfluenzavirus', 'Influenza A virus', 'Influenza A virus (A/Puerto Rico/8/1934(H1N1))'], 'tax_id': 211044, 'title': 'Influenza A virus (A/Puerto Rico/8/1934(H1N1))'},
                      {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Polyomaviridae', 'Betapolyomavirus', 'Human polyomavirus 2', 'JC polyomavirus'], 'tax_id': 10632, 'title': 'JC polyomavirus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Ellioviricetes', 'Bunyavirales', 'Arenaviridae', 'Mammarenavirus', 'Lassa mammarenavirus', 'Lassa mammarenavirus'], 'tax_id': 11620, 'title': 'Lassa mammarenavirus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Rhabdoviridae', 'Rhabdoviridae_u_g', 'Long Island tick rhabdovirus', 'Long Island tick rhabdovirus'], 'tax_id': 1459044, 'title': 'Long Island tick rhabdovirus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Ellioviricetes', 'Bunyavirales', 'Arenaviridae', 'Mammarenavirus', 'Lymphocytic choriomeningitis mammarenavirus', 'Lymphocytic choriomeningitis mammarenavirus'], 'tax_id': 11623, 'title': 'Lymphocytic choriomeningitis mammarenavirus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Filoviridae', 'Marburgvirus', 'Marburg marburgvirus', 'Marburg virus - Musoke, Kenya, 1980'], 'tax_id': 33727, 'title': 'Marburg virus - Musoke, Kenya, 1980'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Paramyxoviridae', 'Morbillivirus', 'Measles morbillivirus', 'Measles morbillivirus'], 'tax_id': 11234, 'title': 'Measles morbillivirus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Flaviviridae', 'Flaviviridae_u_g', 'Menghai flavivirus', 'Menghai flavivirus'], 'tax_id': 1963247, 'title': 'Menghai flavivirus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Poxviridae', 'Molluscipoxvirus', 'Molluscum contagiosum virus', 'Molluscum contagiosum virus subtype 1'], 'tax_id': 10280, 'title': 'Molluscum contagiosum virus subtype 1'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Poxviridae', 'Orthopoxvirus', 'Monkeypox virus', 'Monkeypox virus Zaire-96-I-16'], 'tax_id': 619591, 'title': 'Monkeypox virus Zaire-96-I-16'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Paramyxoviridae', 'Rubulavirus', 'Mumps rubulavirus', 'Mumps rubulavirus'], 'tax_id': 1979165, 'title': 'Mumps rubulavirus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Reoviridae', 'Orbivirus', 'Orbivirus SX-2017a', 'Orbivirus SX-2017a'], 'tax_id': 1955493, 'title': 'Orbivirus SX-2017a'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Poxviridae', 'Parapoxvirus', 'Orf virus', 'Orf virus'], 'tax_id': 10258, 'title': 'Orf virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Picornavirales', 'Picornaviridae', 'Picornaviridae_u_g', 'Pigeon picornavirus B', 'Pigeon picornavirus B'], 'tax_id': 928289, 'title': 'Pigeon picornavirus B'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Rhabdoviridae', 'Lyssavirus', 'Rabies lyssavirus', 'Rabies lyssavirus_u_t'], 'tax_id': 11292, 'title': 'Rabies lyssavirus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Monjiviricetes', 'Mononegavirales', 'Pneumoviridae', 'Pneumoviridae_u_g', 'Respiratory syncytial virus', 'Respiratory syncytial virus'], 'tax_id': 12814, 'title': 'Respiratory syncytial virus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Ellioviricetes', 'Bunyavirales', 'Phenuiviridae', 'Phlebovirus', 'Rift Valley fever phlebovirus', 'Rift Valley fever virus'], 'tax_id': 11588, 'title': 'Rift Valley fever virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Togaviridae', 'Rubivirus', 'Rubella virus', 'Rubella virus'], 'tax_id': 11041, 'title': 'Rubella virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Flaviviridae', 'Flavivirus', 'Saint Louis encephalitis virus', 'Saint Louis encephalitis virus'], 'tax_id': 11080, 'title': 'Saint Louis encephalitis virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Flaviviridae', 'Flavivirus', 'Tick-borne encephalitis virus', 'Tick-borne encephalitis virus_u_t'], 'tax_id': 11084, 'title': 'Tick-borne encephalitis virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Poxviridae', 'Orthopoxvirus', 'Variola virus', 'Variola virus'], 'tax_id': 10255, 'title': 'Variola virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Togaviridae', 'Alphavirus', 'Venezuelan equine encephalitis virus', 'Venezuelan equine encephalitis virus_u_t'], 'tax_id': 11036, 'title': 'Venezuelan equine encephalitis virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Viruses_u_f', 'Viruses_u_g', 'Wenling toga-like virus', 'Wenling toga-like virus'], 'tax_id': 1923542, 'title': 'Wenling toga-like virus'}, {'taxonomy': ['Viruses', 'Viruses_u_p', 'Viruses_u_c', 'Viruses_u_o', 'Togaviridae', 'Alphavirus', 'Western equine encephalitis virus', 'Western equine encephalitis virus'], 'tax_id': 11039, 'title': 'Western equine encephalitis virus'}, {'taxonomy': ['Viruses', 'Negarnaviricota', 'Ellioviricetes', 'Bunyavirales', 'Arenaviridae', 'Mammarenavirus', 'Whitewater Arroyo mammarenavirus', 'Whitewater Arroyo mammarenavirus'], 'tax_id': 46919, 'title': 'Whitewater Arroyo mammarenavirus'}],
sample_metadata=[{'name': 'GCF_000861845.1_ViralProj15432_genomic.fna', 'reads_total': 1, 'label': 'label0', 'label_name': 'No Label', 'file_id': '4cc50293-3ac9-4565-98d2-0395f5ff14d1', 'dataset_id': u'e0d0e947-1d20-4ee0-918c-c763c429f8ad'}, {'name': u'GCF_000861845.1_ViralProj15432_genomic.fna', 'reads_total': 1, 'label': 'label0', 'label_name': 'No Label', 'file_id': '29c6d2e6-617c-494c-b868-6b6ec5518d81', 'dataset_id': '9c27b10c-687e-43d2-9631-21a416489541'}, {'name': 'virus_setB_10K.fastq.fasta', 'reads_total': 390003, 'label': 'label0', 'label_name': 'No Label', 'file_id': '1af0ba9b-b94a-485e-9232-37dffc48754f', 'dataset_id': 'db1518c1-8779-423b-9f48-92d1c2bb77a3'}],
table_id='05775205-8479-4346-9866-a45bdb449d70',
type = "OTU table")

#
# output = (array([[ 1.        ,  0.        , -0.        ],
#                  [42.        ,  0.97297728,  5.24199479]]),
#           array([u'GCF_000861845.1_ViralProj15432_genomic.fna', u'virus_setB_10K.fastq.fasta'], dtype=object),
#           ['chao1', 'simpson', 'shannon'])

data_3 = array([[0.,     0.,        5.29367628],
 [0.,         0.,         5.52763892],
 [0.,        0.,         1.43136376],
 [0.,         0.,         5.55793407],
 [0. ,        0.  ,       5.53187088],
 [0.  ,       0.  ,       1.69897   ],
 [0.   ,      0.     ,    5.56313327],
 [0.    ,     0.     ,    5.11045839],
 [3.53096768 , 3.53096768 , 0.],
 [0. ,        0. ,        5.48072107],
 [0.  ,       0.   ,      5.45509927],
 [0.     ,    0.  ,       5.50020628],
 [0.  ,       0. ,        5.50599603],
 [0.     ,    0.     ,    5.59733532],
 [0. ,        0.  ,       1.11394335],
 [0.  ,       0.   ,      5.58181886],
 [0.  ,       0. ,        5.51705331],
 [0.  ,       0. ,        5.42805538],
 [0.  ,       0.  ,       5.49400981],
 [0.   ,      0.   ,      5.56329831],
 [0.   ,      0.  ,       5.41639084],
 [0. ,        0.   ,      5.5121452 ],
 [0.  ,       0.  ,       5.56240945],
 [0.   ,      0.  ,       5.58252697],
 [0.   ,      0.  ,       5.52382058],
 [0.  ,       0.      ,   5.50203362],
 [0.  ,       0.    ,     5.39554654],
 [0.  ,       0.  ,       5.1721766 ],
 [0.  ,       0.     ,    5.55956312],
 [0.     ,    0.  ,       5.42083006],
 [0.   ,      0.  ,       5.29768831],
 [0.   ,      0.  ,       5.53319355],
 [0.   ,      0.    ,     5.54272586],
 [0.  ,       0.     ,    5.61124096],
 [0.    ,     0.  ,       5.54092606],
 [0.    ,     0.   ,      5.31292593],
 [0.  ,       0.    ,     5.50991456],
 [0.    ,     0.  ,       5.44642842],
 [0.  ,       0.    ,     5.12015589],
 [0.    ,     0.  ,       5.50973867],
 [0.  ,       0.       ,  5.4030879 ],
 [0.  ,       0.  ,       5.51752842],
 [0.   ,      0.   ,      5.57927586]])





# print(data.transpose()[2])
# diver = chao1(data.transpose()[2])
# diver_shannon = shannon(data.transpose()[2])
# diver_simpson = simpson(data.transpose()[2])
# print(diver)
# print(diver_shannon)
# print(diver_simpson)
#

# biom = getBiomData(input_biom)
# print(type(biom.data))
ALPHA_DIVERSITY_METHODS = ['chao1', 'simpson', 'shannon']


def alpha(table: biom.Table):
    if table.is_empty():
        raise ValueError("The provided table object is empty")
    alpha_diversities = []
    for metric in ALPHA_DIVERSITY_METHODS:
        counts = table.matrix_data.toarray().astype(float).T
        sample_ids = table.ids(axis='sample')

        result = skbio.diversity.alpha_diversity(metric=metric, counts=counts,
                                                 ids=sample_ids)
        result.name = metric
        print(result)
        print(type(result))
        alpha_diversities.append(result)
    return alpha_diversities


def _format_alpha_results(result, sample_metadata=None):
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

res = alpha(input_biom)
# r1 = _format_alpha_results(input_biom, sample_metadata=)

# chao1_div = alpha_diversity('chao1', data.transpose()[2])
# shannon_div = alpha_diversity('shannon', data.transpose()[2])
# simpson_div = alpha_diversity('simpson', data.transpose()[2])
# #
# print(colored('\n Chao1 diversity: \n', 'green'))
# print(chao1_div)
# print(colored('\n Shannon diversity: \n', 'green'))
# print(shannon_div)
# print(colored('\n Simpson diversity: \n', 'green'))
# print(simpson_div)

