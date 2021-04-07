#!/usr/bin/env python

from rdkit import Chem
import os
import csv
import numpy as np


def csv_read():
    """
    读csv文件，提取指定列
    """
    file_location = 'D:/test/GLORY_test_dataset.csv'

    with open(file_location, 'r', encoding='utf-8') as csv_file:

        SUB_SMILES_name = 0
        PRO_SMILES_name = 0
        for row in csv_file:
            SUB_SMILES = str(row['parent SMILES'])
            PRO_SMILES = str(row['PRO_SMILES'])
            SUB_SMILES_name += 1
            PRO_SMILES_name += 1

            if 'None' in SUB_SMILES:
                with open('D:/test/GLORY_SUB_mol/' + "{}".format(SUB_SMILES_name) + 'None.mol', 'a') as txt_file:
                    txt_file.write('None')
            else:
                SUB_mol = Chem.MolFromSmiles(SUB_SMILES)
                if 'None' in str(SUB_mol):
                    with open('D:/test/GLORY_SUB_mol/'+"{}".format(SUB_SMILES_name) + 'None.mol', 'a') as txt_file:
                        txt_file.write('None')
                else:
                    Chem.MolToMolFile(SUB_mol, 'D:/test/GLORY_SUB_mol/' + "{}".format(SUB_SMILES_name) + '.mol')

            if 'None' in PRO_SMILES:
                with open('D:/test/Drugbank_PRO_mol/' + "{}".format(PRO_SMILES_name) + 'None.mol', 'a') as txt_file:
                    txt_file.write('None')
            else:
                PRO_mol = Chem.MolFromSmiles(PRO_SMILES)
                if 'None' in str(PRO_mol):
                    with open('D:/test/Drugbank_PRO_mol/'+"{}".format(PRO_SMILES_name)+ 'None.mol', 'a') as txtt_file:
                        txtt_file.write('None')
                else:
                    Chem.MolToMolFile(PRO_mol, 'D:/test/Drugbank_PRO_mol/' + "{}".format(PRO_SMILES_name) + '.mol')

            with open('D:/test/Drugbank_SUB_SMILES/'+"{}".format(SUB_SMILES_name)+ '.txt', 'a') as txt_file:
                txt_file.write(SUB_SMILES)
            txt_file.close()
            with open('D:/test/Drugbank_PRO_SMILES/'+"{}".format(SUB_SMILES_name) + '.txt', 'a') as PRO_file:
                PRO_file.write(PRO_SMILES)
            PRO_file.close()

class CsvToRdt(object):

    """
    csv to rdt_input
    """
    def __init__(self, input_file, output_file):

        self.input_file_location = input_file
        self.output_file_location = output_file

    def read(self):

        with open(self.input_file_location, 'r', encoding='utf-8') as csv_file:
            csv_file = csv.DictReader(csv_file)
            n = 0
            with open(self.output_file_location, 'a') as txt_file:
                for row in csv_file:
                    SUB_SMILES = str(row['parent SMILES'])
                    PRO_SMILES = str(row['metabolite SMILES (separated by space)']).split(' "')
                    for i in range(len(PRO_SMILES)):
                        txt_file.write(SUB_SMILES + '>>' + PRO_SMILES[i].strip('"') + '\n')
                        n += 1
            txt_file.close()
            print(n)

    def fame3(self):
        writer = Chem.SDWriter(self.output_file_location)
        writer.SetProps(['Name'])
        with open(self.input_file_location, 'r', encoding='utf-8') as csv_file:
            csv_file = csv.DictReader(csv_file)
            mols = []
            name = []
            for row in csv_file:
                PRO_SMILES = str(row['metabolite SMILES (separated by space)']).split(' "')
                for j in range(len(PRO_SMILES)):
                    name.append(str(row['nameofparentmolecule']))
                    SUB_SMILES = str(row['parent SMILES'])
                    mols.append(Chem.MolFromSmiles(SUB_SMILES))
            for i, mol in enumerate(mols):
                mol.SetProp('Name', name[i])
                writer.write(mol)
        writer.close()


if __name__ == '__main__':

    # csv_location =  ('D:/test/GLORY_test_dataset.csv')
    # txt_location = ('D:/test/GLORY_rdt_input/GLORY_rdt_input.txt')
    ctr = CsvToRdt('D:/test/GLORY_test_dataset.csv', 'D:/test/GLORY_rdt_input/GlORY.sdf')
    csv_read()
    # ctr.fame3()
