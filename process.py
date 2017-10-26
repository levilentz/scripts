#!/global/homes/l/llentz/anaconda3/bin/python

import tempfile
import sys
import glob
import QE
from bandgapoccu import bg as BG
import os
import json
import re

def separate_aimd(fname='scf.out',tempfile='OUT.out',jfile='data.json'):
  temp = [zz.split()[2] for zz in open(fname).read().split('\n') if 'temperature' in zz]
  temp = [float(zz) for zz in temp if re.match('\d+.\d+$',zz)]
  scf = open(fname,'r')

  header = ''

  for i in scf:
    if 'PseudoPot' in i: break
    header += i


  data = []
  rest = header + '\n\n'
  tmp = ''
  tmp += rest
  cnt = 1
  try:
    for i in scf:
      while '!    total energy' not in i:
        i = next(scf)
        tmp += i
      #data[cnt] = {}
      t_ = {}
      OUT = tempfile
      output = open(OUT,'w')
      output.write(tmp)
      output.close()
      x = QE.Struct()
      x.File_Process(OUT)
      atoms,cell = x.return_params()
      t_['str'] = x.print()
      t_['atoms'] = atoms
      t_['cell'] = cell
      del x
      Gap = BG(OUT)
      t_['gap'] = Gap.bg
      try:
        t_['temp'] = temp[cnt - 1] #[zz for zz in open(str(cnt) + '.out').read().split('\n') if 'temperature' in zz][-1].split()[2]
      except:
        t_['temp'] = 'NA'
      cnt += 1
      tmp = rest
      data.append(t_)
      with open(jfile, 'w') as outfile:
          json.dump(data, outfile)
      i = next(scf)
      continue
    scf.close()
  except StopIteration:
    scf.close()
    with open(jfile, 'w') as outfile:
        json.dump(data, outfile)

if __name__ == '__main__':
  separate_aimd(fname=sys.argv[1])
