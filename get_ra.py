from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import pandas as pd
import numpy as np
import os
from tqdm import tqdm

class GetReaAtom(object):
    def __init__(self,aam_file, w):
        self.aam_file = aam_file
        self.w = w

    @staticmethod
    def map2id(mol):
        dic = {}
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map == 0:
                pass
                # print('atom {}'.format(atom.GetSymbol()), 'id: {} not map'.format(atom.GetIdx()))
            else:
                dic[atom.GetAtomMapNum()] = atom.GetIdx()
        return dic

    @staticmethod
    def id2map(mol):
        dic = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                dic[atom.GetIdx()] = atom.GetAtomMapNum()
        return dic

    def sub_map2pro_map2pro_id(self, sub, pro):
        """input:sub_mol, pro_mol  return:sub_pro_id"""
        dic = {}
        promap2id = self.map2id(pro)
        for sub_atom in sub.GetAtoms():
            sub_map = sub_atom.GetAtomMapNum()  # sub map
            if sub_map in promap2id.keys():  # 确保原子被map
                dic[sub_map] = promap2id[sub_map]
        return dic

    def pro_map2sub_map2sub_id(self, sub, pro):
        """input:sub_mol, pro_mol  return:pro_sub_id"""
        dic = {}
        submap2id = self.map2id(sub)
        for pro_atom in pro.GetAtoms():
            pro_map = pro_atom.GetAtomMapNum()  # pro map
            if pro_map in submap2id.keys():
                dic[pro_map] = submap2id[pro_map]
        return dic

    @staticmethod
    def compare_adj(react_atom, sub, pro, subid2map, submap2id, promap2id, proid2map, sub_pro_id, pro_sub_id):
        for i in range(sub.shape[0]):
            for j in range(sub.shape[1]):
                s1 = sub[i, j]  # sub adj的值
                if i in subid2map.keys() and j in subid2map.keys():  # 不考虑没有map的原子
                    if_i = subid2map[i] in promap2id.keys()  # 判断 i和j对应的map是否在产物里
                    if_j = subid2map[j] in promap2id.keys()
                    if if_i == 1 and if_j == 1:  # 两端都在生成物里
                        p1 = pro[sub_pro_id[subid2map[i]], sub_pro_id[subid2map[j]]]
                        if s1 != 0 and s1 != p1:  # 键的类型改变/键断开,而且两端都在生成物里,则i和j是反应位点
                            react_atom.add(i)
                            react_atom.add(j)
                    elif if_i + if_j == 1 and s1 != 0:  # 键断开,有一端不在产物里,则i和j是反应位点
                        react_atom.add(i)
                        react_atom.add(j)
        if len(react_atom) == 0:
            for m in range(pro.shape[0]):
                for n in range(pro.shape[1]):
                    p2 = pro[m, n]  # sub adj的值
                    if m in proid2map.keys() and n in proid2map.keys() and p2 != 0:  # 不考虑没有map的原子
                        if_m = proid2map[m] in submap2id.keys()  # 判断 i和j对应的map是否在反应物里
                        if_n = proid2map[n] in submap2id.keys()
                        if if_m == 1 and if_n == 1:  # 两端都在反应物物里
                            s2 = sub[pro_sub_id[proid2map[m]], pro_sub_id[proid2map[n]]]
                            if  s2 == 0:  # 新键生成,则i和j对应的sub的id是反应位点
                                react_atom.add(pro_sub_id[proid2map[m]])
                                react_atom.add(pro_sub_id[proid2map[n]])
                        elif if_m == 0 and if_n == 1:  # 键生成,有一端不在产物里,则在反应物里的对应的sub的id是反应位点
                            react_atom.add(pro_sub_id[proid2map[n]])
                        elif if_m == 1 and if_n == 0:
                            react_atom.add(pro_sub_id[proid2map[m]])
        # print('{}是反应位点'.format(reaction_atom))
        return react_atom

    def get_mol(self):
        rxn = open(self.aam_file, 'r').readlines()
        sub_mol = []
        pro_mol = []
        nameall = []
        for line in rxn:
            sub_smi = line.split('>>')[1]
            pro_smi = line.split('>>')[2]
            name = line.split('_ECBLAST')[0]
            sub_mol.append(Chem.MolFromSmiles(sub_smi))
            pro_mol.append(Chem.MolFromSmiles(pro_smi))
            nameall.append(name)
        return sub_mol, pro_mol, nameall

    def get_all(self):
        sub_mol, pro_mol, name = self.get_mol()
        for i in tqdm(range(len(sub_mol))):
            react_atom = set()
            sub_adj = rdmolops.GetAdjacencyMatrix(sub_mol[i], useBO=1).astype(int)
            pro_adj = rdmolops.GetAdjacencyMatrix(pro_mol[i], useBO=1).astype(int)
            sub_map2id = self.map2id(sub_mol[i])
            pro_map2id = self.map2id(pro_mol[i])
            # print('sub_id', sub_map2id)
            # print('pro_id', pro_map2id)
            sub_id2map = self.id2map(sub_mol[i])
            pro_id2map = self.id2map(pro_mol[i])
            # print('sub_map', sub_id2map)
            # print('pro_map', pro_id2map)
            sub2pro_id = self.sub_map2pro_map2pro_id(sub_mol[i], pro_mol[i])
            pro2sub_id = self.pro_map2sub_map2sub_id(sub_mol[i], pro_mol[i])
            # print('sub_pro_id', sub2pro_id)
            # print("name", name[i])
            ra = self.compare_adj(react_atom, sub_adj, pro_adj, sub_id2map, sub_map2id, pro_map2id,
                                  pro_id2map, sub2pro_id, pro2sub_id)
            res_str = ''
            if len(ra) > 0:
                for j in ra: res_str = res_str + ' ' + str(j)
                sub_mol[i].SetProp('SOM', res_str)
                sub_mol[i].SetProp('_Name', name[i])
                w.write(sub_mol[i])
            # print('--------------------------------------------------------')
        print('done')


if __name__ == "__main__":
    aam = os.path.join(r'D:\aam_data\data\Met\aam\mdb_aam.txt')
    w = Chem.SDWriter(r'D:\aam_data\data\Met\mdb.sdf')
    GRA = GetReaAtom(aam, w)
    GRA.get_all()