#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 14:07:26 2023

@author: maria
"""
## paper 4
# change name genes

import os

path = os.path.dirname(os.path.realpath('changeNameGene_paper4.py'))

pruebas = []

f = open(path + '/normalBreast_DEGs_table_symbol.csv','w')
degs = open(path + '/normalBreast_DEGs_table.csv','r')
for line in degs:
    if line.startswith('attr_name'):
        f.write(line)
    else:
        line = line.strip()
        spl = line.split(',')
        activated = open(path + '/activatedGenes_breast_table.txt','r')
        for line1 in activated:
            if spl[0] in line1:
                spl1 = line1.split('\t')
                if not spl1[1] in pruebas:
                    pruebas.append(spl1[1])
                    f.write(spl1[1] + "," + ",".join(spl[1:]) + "\n")
                else:
                    pass
        activated.close()
        repressed = open(path + '/repressedGenes_breast_table.txt','r')
        for line1 in repressed:
            if spl[0] in line1:
                spl1 = line1.split('\t')
                if not spl1[1] in pruebas:
                    pruebas.append(spl1[1])
                    f.write(spl1[1] + "," + ",".join(spl[1:]) + "\n")
                else:
                    pass
        repressed.close()
degs.close()
f.close()

f = open(path + '/tumorBreast_DEGs_table_symbol.csv','w')
degs = open(path + '/tumorBreast_DEGs_table.csv','r')
for line in degs:
    if line.startswith('attr_name'):
        f.write(line)
    else:
        line = line.strip()
        spl = line.split(',')
        activated = open(path + '/activatedGenes_breast_table.txt','r')
        for line1 in activated:
            if spl[0] in line1:
                spl1 = line1.split('\t')
                if not spl1[1] in pruebas:
                    pruebas.append(spl1[1])
                    f.write(spl1[1] + "," + ",".join(spl[1:]) + "\n")
                else:
                    pass
        activated.close()
        repressed = open(path + '/repressedGenes_breast_table.txt','r')
        for line1 in repressed:
            if spl[0] in line1:
                spl1 = line1.split('\t')
                if not spl1[1] in pruebas:
                    pruebas.append(spl1[1])
                    f.write(spl1[1] + "," + ",".join(spl[1:]) + "\n")
                else:
                    pass
        repressed.close()
degs.close()
f.close()
