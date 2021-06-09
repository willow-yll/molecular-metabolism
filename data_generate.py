from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import os
import numpy as np
import tqdm


"""
data for GCC
https://github.com/THUDM/GCC
"""
file_dir_path = r'D:\test_data'
file_name = r'glory_aam.sdf'
file_path = os.path.join(file_dir_path,file_name)
suppls = Chem.SDMolSupplier(file_path)


def name():
    f_name = open(os.path.join(r'D:\test_data\data', 'name.txt'), 'w')
    for mol in suppls:
        name = mol.GetProp('_Name')
        f_name.writelines(name+'\n')
    f_name.close()


def node_label():
    for mol in suppls:
        name = mol.GetProp('_Name')
        f = open(os.path.join(r'D:\test_data\data','{}.nodelabel'.format(name)), 'w')
        som = mol.GetProp('SOM').split()
        nodes = [a.GetIdx() for a in mol.GetAtoms()]
        label = [0 for i in range(len(nodes))]
        for i in som:
            label[int(i)] = 1
        for j in range(len(nodes)):
            f.writelines(str(nodes[j])+' '+str(label[j])+'\n')
        f.close()


def edge_list():
    for mol in suppls:
        name = mol.GetProp('_Name')
        f = open(os.path.join(r'D:\test_data\data', '{}.edgelist'.format(name)), 'w')
        adj = Chem.rdmolops.GetAdjacencyMatrix(mol)
        a = []
        b = []
        for i in range(len(adj)):
            for j in range(len(adj)):
                if adj[i][j] == 1:
                    a.append(i)
                    b.append(j)
        for k in range(len(a)):
            f.writelines(str(a[k]) + ' ' + str(b[k]) + '\n')
        f.close()


def all_node_label():
    nodes = []
    labels = []
    f_n = open(os.path.join(r'D:\test_data\data', 'dm_glory.nodelabel'), 'w')
    f_e = open(os.path.join(r'D:\test_data\data', 'dm_glory.edgelist'), 'w')
    a_s = []
    b_s = []
    for mol in suppls:
        som = mol.GetProp('SOM').split()
        adj = Chem.rdmolops.GetAdjacencyMatrix(mol)
        a = []
        b = []
        for i in range(len(adj)):
            for j in range(len(adj)):
                if adj[i][j] == 1:
                    a.append(i)
                    b.append(j)
        a = [m + len(nodes) for m in a]
        b = [n + len(nodes) for n in b]
        a_s.extend(a)
        b_s.extend(b)
        node = [a.GetIdx() for a in mol.GetAtoms()]
        node = [b+len(nodes) for b in node]
        nodes.extend(node)
        label= [0 for i in range(len(node))]
        for c in som:
            label[int(c)] = 1
        labels.extend(label)
    for d in range(len(nodes)):
        f_n.writelines(str(nodes[d])+' '+str(labels[d])+'\n')
    f_n.close()
    for e in range(len(a_s)):
        f_e.writelines(str(a_s[e]) + ' ' + str(b_s[e]) + '\n')
    f_e.close()



if __name__ == "__main__":
    all_node_label()
