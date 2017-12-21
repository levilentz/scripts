#!/usr/bin/env python

import QE
from lxml import etree as ET
from random import shuffle
import sys


def xml(f,extra_tag=None):
  '''This takes a list of directories and creates a PROPhet xml file'''
  root = ET.Element("PROPhet")
  nsystms = ET.Element('nsystem')
  nsystms.text = str(len(f))
  root.append(nsystms)
  systems = ET.Element("systems")
  for i in range(200): shuffle(f)
  N_train = int(0.80*len(f))
  N_val = int(0.90*len(f))
  sys = []
  for c,i in enumerate(f):
    if c%10 == 0 : print(c)
    if c < N_train: t_flag = "train"
    elif N_train < c < N_val: t_flag = "val"
    else: t_flag = 'test'
    x = QE.Struct()
    x.XML_Process(i)
    system = ET.Element("system",id=str(c + 1))
    train = ET.Element('train')
    train.text = t_flag
    system.append(train)
    lattice = ET.Element('lattice',units='angstrom')
    for j in x.lattice:
      l = ET.Element(j)
      l.text = ' '.join([str(zz) for zz in x.lattice[j]])
      lattice.append(l)
    system.append(lattice)
    atoms = ET.Element('atoms',units='angstrom')
    atm,cell = x.return_params()
    natoms = ET.Element('natoms')
    natoms.text = str(len(x.atoms))
    species = ET.Element('species')
    ntype = len(set([zz.split()[0] for zz in atm.split('\n')[:-1]]))
    species.text = str(ntype)
    atoms.append(natoms)
    atoms.append(species)
    for i in atm.split('\n')[:-1]: 
      atom = ET.Element("atom",specie=i.split()[0])
      atom.text = ' '.join(i.split()[1:4])
      atoms.append(atom)
    system.append(atoms)
    target = ET.Element('target')
    target.text = str(x.energy)
    system.append(target)
    if extra_tag is not None:
      tag = extra_tag[c]['tag']
      val = extra_tag[c]['val']
      if 'other_tags' in list(extra_tag[c].keys()):
        _ = Et.Element(tag,**extra_tag[c])
      else: 
        _ = Et.Element(tag)
      _.text = val
      system.append(_)
    sys.append(system)
    del x
  for i in sys:
    systems.append(i)
  root.append(systems)
  str_ = ET.tostring(root,pretty_print=True).decode('utf-8')
  return str_

if __name__ == '__main__':
  d = open(sys.argv[1]).read().split()[0:10]
  t = xml(d)
  f = open('PROPhet.xml','w')
  f.write(t)
  f.close()
