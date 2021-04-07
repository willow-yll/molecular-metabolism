# -*- coding: utf-8 -*-
# @Time     : 2021/3/15 17:56
# Author    : YangLiu
# @Site     : 
# @File     : read_json.py

import demjson
import json
from rdkit import Chem
import os

class JsonRead:

    """
    把smarts写成mol
    """
    @staticmethod
    def read(file_path, new_file):
        with open(file_path, 'r', encoding='utf8') as f:
            json_to_dict = json.load(f) #list 反应对
            n = 0
            for i in range(len(json_to_dict)):
                pro_dic = json_to_dict[i]['Metabolites']
                n += 1
                if type(pro_dic) == dict:    #dict
                    pro = pro_dic['SMILES']
                    sub_dic = json_to_dict[i]['Parent molecule']  #dict
                    sub = sub_dic['SMILES']
                    with open(os.path.join(new_file, "{}".format(n) , '.txt'), 'a') as txt_file:
                        txt_file.write(sub+'>>'+pro)
                    txt_file.close()
                else:
                    m = 0
                    for j in pro_dic: #j 是list中的dict
                        m += 1
                        pro = j['SMILES']
                        sub_dic = json_to_dict[i]['Parent molecule']  # dict
                        sub = sub_dic['SMILES']
                        with open(os.path.join(new_file, "{}".format(n)+'_'+"{}".format(m)+'.txt'), 'a') as txt_file:
                            txt_file.write(sub + '>>' + pro)
                        txt_file.close()
        f.close()

    @staticmethod
    def inchi(file_path, new_file):
        f = open(file_path, 'r', encoding='utf8')
        json_to_dict = json.load(f)
        smarts = json_to_dict['biotransformations']
        for key,value in smarts.items():
            n = 0
            d = 0
            pro_list = value['Products']
            for i in pro_list:
                d += 1
                pro_InChI = i['InChI']
                pro_m = Chem.inchi.MolFromInchi(pro_InChI)
                sub_InChI = value['Substrate']['InChI']
                m = Chem.inchi.MolFromInchi(sub_InChI)
                try:
                    if d>1 :
                        n += 1
                        wf = open(os.path.join(new_file, key.split('D')[1]+'_'+'{}'.format(n)+'.txt'), 'w', encoding='utf8')
                        sub_smi = Chem.MolToSmiles(m)
                        pro_smi = Chem.MolToSmiles(pro_m)
                        wf.write(key.split('D')[1]+'>>'+sub_smi+'>>'+pro_smi)
                        wf.close()
                    else:
                        wf = open(os.path.join(new_file, key.split('D')[1]+'.txt'), 'w',
                                  encoding='utf8')
                        sub_smi = Chem.MolToSmiles(m)
                        pro_smi = Chem.MolToSmiles(pro_m)
                        wf.write(key.split('D')[1] + '>>' + sub_smi + '>>' + pro_smi)
                        wf.close()
                except Exception as e:
                    print(key,e)
        f.close()

class SmirkToRdt:

    """
    把smirks写成rdt输入
    """

    @staticmethod
    def read(json_path, out_path):

        json_file = open(json_path, 'r', encoding='utf8').read()
        dict_file = demjson.decode(json_file)
        reactions = dict_file['reactions']

        with open(out_path, 'a') as txt_file:
            n = 0
            for key, value in reactions.items():
                smirks_str = "".join(value['smirks'])
                if len(smirks_str.split(">>", 1)) == 2:
                    n += 1
                    sub = smirks_str.split(">>", 1)[0]
                    sub = Chem.MolFromSmarts(sub)
                    sub = Chem.MolToSmiles(sub)
                    pro = smirks_str.split(">>", 1)[1]
                    pro = Chem.MolFromSmarts(pro)
                    pro = Chem.MolToSmiles(pro)
                    txt_file.write(sub + '>>' + pro + '\n')
                else:
                    txt_file.write('\n')
            print(n)
        txt_file.close()

    @staticmethod
    def fame3(input_path, output_path):
        pass


if __name__ == "__main__":

    # json_path = 'D:/test/MetabolicReactions.json'
    # out_path = 'D:/test/MetabolicReactions_rdt_input/MetabolicReactions_rdt_input_smi.txt'
    # s_to_r = SmirkToRdt()
    # s_to_r.read(json_path, out_path)
    # filepath = 'D:/test/MetabolicReactions.json'
    # filepath = "/data/yangliu/dataset/GLORY_reference_dataset.json"
    # new = "/data/yangliu/dataset"
    metdata_new = "/data/yangliu/dataset/MetXBioDB"
    METDB = "/data/yangliu/dataset/MetXBioDB/MetXBioDB-1-0.json"
    jr = JsonRead()
    jr.inchi(METDB, metdata_new)
    # jr.read(filepath, new)