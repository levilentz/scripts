#!/usr/bin/env python

import tempfile
import sys
import glob
import QE
from bandgapoccu import bg as BG
import os
import json
import io
import re
temp = [zz.split()[2] for zz in open(sys.argv[1]).read().split('\n') if 'temperature' in zz]
temp = [float(zz) for zz in temp if re.match('\d+.\d+$',zz)]
scf = open(sys.argv[1],'r')

header = ''

for i in scf:
  if 'PseudoPot' in i: break
  header += i


data = {}
rest = header + '\n\n'
tmp = ''
tmp += rest
cnt = 1
try:
  for i in scf:
    while '!    total energy' not in i:
      i = next(scf)
      tmp += i
    #if not os.path.exists('%s.save'%cnt): os.makedirs('%s.save'%cnt)
    data[cnt] = {}
    OUT = 'OUT.out'
    output = open(OUT,'w')
    output.write(tmp)
    output.close()
    x = QE.Struct()
    x.File_Process(OUT)
    atoms,cell = x.return_parms()
    data[cnt]['str'] = x.print()
    data[cnt]['atoms'] = atoms
    data[cnt]['cell'] = cell
    del x
    Gap = BG(OUT)
    data[cnt]['gap'] = Gap.bg
    try:
      data[cnt]['temp'] = temp[cnt - 1] #[zz for zz in open(str(cnt) + '.out').read().split('\n') if 'temperature' in zz][-1].split()[2]
    except:
      data[cnt]['temp'] = 'NA'
    cnt += 1
    tmp = rest
    with open('data.json', 'w') as outfile:
        json.dump(data, outfile)
    i = next(scf)
    continue
  scf.close()
except:
  scf.close()
  with open('data.json', 'w') as outfile:
      json.dump(data, outfile)
