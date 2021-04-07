# -*- coding: utf-8 -*-
# @Time     : 2021/3/12 13:48
# Author    : YangLiu
# @Site     : 
# @File     : Draw.py
# @Software :


from rdkit.Chem import Draw,AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
# from rdkit.Chem.Draw import IPythonConsole
from rdkit import Chem
import os
# IPythonConsole.ipython_useSVG = True


class DrawMolWithIdx:

    """
    从mol文件绘制带标签的原子图
    """

    def __init__(self):
        pass

    @staticmethod
    def mol_with_atom_index(mol_file):

        atoms = mol_file.GetNumAtoms()
        for idx in range(atoms):
            mol_file.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))

    @staticmethod
    def draw_mol(mol_file, save_mol_path):

        opts = DrawingOptions()
        opts.includeAtomNumbers = True
        opts.includeAtomNumbers = True
        opts.bondLineWidth = 2.8
        Draw.MolToFile(
            mol_file,  # mol对象
            save_mol_path,  # 图片存储地址
            size=(4000, 2000),
            kekulize=True,
            wedgeBonds=True,
            imageType=None,
            fitImage=False,
            options=opts,
        )


class DrawReaction:
    """
    绘制反应图
    """

    def __init__(self):
        pass

    @staticmethod
    def draw(reaction, save_reaction_path):
        d2d = Draw.MolDraw2DCairo(6000, 2000)
        d2d.DrawReaction(reaction)
        png = d2d.GetDrawingText()
        # 'D:/test/image/output.png'
        open(save_reaction_path, 'wb+').write(png)

    @staticmethod
    def draw_highlight(reaction, save_reaction_path):
        """
        high light版
        """
        d2d = Draw.MolDraw2DCairo(4000, 2000)
        d2d.DrawReaction(reaction, highlightByReactant=True)
        png = d2d.GetDrawingText()
        open(save_reaction_path, 'wb+').write(png)

if __name__ == '__main__':

    # """
    # 绘制mol图
    # """
    # # mol = Chem.MolFromSmarts('COC(=O)[C@@H]1CC2=CC(=O)CC[C@]2(C)[C@@]23O[C@@H]2C[C@@]2(C)[C@@H](CC[C@@]22CC(=O)C(=O)O2)[C]13COC(=O)[C@@H]1CC2=CC(=O)CC[C@]2(C)[C@@]23O[C@@H]2C[C@@]2(C)[C@@H](CC[C@@]22CC(=O)C(=O)O2)[C@H]13')
    # mol = Chem.MolFromSmarts('[*:1][C:2]=1[C:3]([H:9])=[C:5]([H:11])[C:6](=[C:4]([H:10])[C:7]1[H:12])[*:8]')
    # mol_path = 'D:/test/image/test_1.png'
    # DrawIdx = DrawMolWithIdx()
    # DrawIdx.mol_with_atom_index(mol)
    # DrawIdx.draw_mol(mol, mol_path)

    """
    绘制reaction图
    """
    rxn = AllChem.ReactionFromSmarts('[O:1]=[C:2]1[N:3]([N:4]=[C:5]2[CH:6]=[CH:7][CH:8]=[CH:9][N:10]12)[CH2:11][CH2:12][CH2:13][N:14]3[CH2:15][CH2:16][N:17]([C:18]=4[CH:19]=[CH:20][CH:21]=[C:22]([Cl:23][CH:24]5[N:25]([C:26]=6[CH:27]=[CH:28][CH:29]=[C:30]([Cl:31][CH:54]7[N:55]([C:56]=8[CH:57]=[CH:58][CH:59]=[C:60]([Cl:61])[CH:62]8)[CH2:63][CH2:64][N:65]([CH2:66][CH2:67][CH2:68][N:69]9[N:70]=[C:71]%10[CH:72]=[CH:73][CH:74]=[CH:75][N:76]%10[C:77]9=[O:32])[CH2:78]7)[CH:33]6)[CH2:34][CH2:35][N:36]([CH2:37][CH2:38][CH2:39][N:40]%11[N:41]=[C:42]%12[CH:43]=[CH:44][CH:45]=[CH:46][N:47]%12[C:48]%11=[O:49])[CH2:50]5)[CH:51]4)[CH2:52][CH2:53]3>>[O:1]=[C:2]1[N:3]([N:4]=[C:5]2[N:10]1[CH:9]=[CH:8][CH:7]3[O:32][CH:6]23)[CH2:11][CH2:12][CH2:13][N:14]4[CH2:15][CH2:16][N:17]([CH2:52][CH2:53]4)[C:18]=5[CH:51]=[C:22]([CH:21]=[CH:20][CH:19]5)[Cl:23][CH:24]6[N:25]([C:26]7=[CH:27][CH:28]=[CH:29][C:30]([Cl:31])=[CH:33]7)[CH2:34][CH2:35][N:36]([CH2:37][CH2:38][CH2:39][N:40]8[N:41]=[C:42]9[N:47]([CH:46]=[CH:45][CH:44]%10[O:79][CH:43]9%10)[C:48]8=[O:49])[CH2:50]6')
    DrawReaction = DrawReaction()
    png_path = '/data/yangliu/dataset/1012.png'
    DrawReaction.draw(rxn, png_path)




