#!/usr/bin/env python

import QE
from lxml import etree as ET
from random import shuffle
import sys


def xml(f,extra_tag=None,p_=10):
  '''This takes a list of directories and creates a PROPhet xml file'''
  root = ET.Element("PROPhet")
  nsystms = ET.Element('nsystem')
  root.append(nsystms)
  systems = ET.Element("systems")
  for i in range(200): shuffle(f)
  N_train = int(0.80*len(f))
  N_val = int(0.90*len(f))
  sys = []
  cnt = 1
  for c,i in enumerate(f):
    if c%p_ == 0 : print(c)
    if c < N_train: t_flag = "train"
    elif N_train < c < N_val: t_flag = "val"
    else: t_flag = 'test'
    x = QE.Struct()
    try:
      x.XML_Process(i)
    except:
      continue
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
        _ = ET.Element(tag,**extra_tag[c]['other_tags'])
      else: 
        _ = ET.Element(tag)
      _.text = val
      system.append(_)
    sys.append(system)
    cnt += 1
    del x
  for i in sys:
    systems.append(i)
  nsystms.text = str(cnt)
  root.append(systems)
  str_ = ET.tostring(root,pretty_print=True).decode('utf-8')
  return str_

if __name__ == '__main__':
  d = open(sys.argv[1]).read().split()[0:10]
  t = xml(d)
  f = open('PROPhet.xml','w')
  f.write(t)
  f.close()
