import pandas as pd
import json
import numpy as np
import QE 
from bandgapoccu import bg
import json, os

def get_hse(all_,phse,hse_location,pbe_fname='.out',hse_fname='.out'):
  '''This matches the HSE data with the PBE data'''
  t = [i for i in all_ if i['phase'] == phse]
  hse_lst = []
  for c,i in enumerate(t):
    key = i['location'].split('/')[-1].replace('.save','')
    d_ = i['location'][0:i['location'].rfind('/')]
    pbe = bg(d_ + '/' + key + pbe_fname)
    i['bandgap'] = pbe.bg
    i['metallic'] = pbe.metallic
    hse_f = '/'.join([hse_location,i['dopant'],key+hse_fname])
    hse = None
    if os.path.isfile(hse_f):
      if 'job done' in open(hse_f).read().lower():
        hse = bg(hse_f)
        i['hse_bandgap'] = hse.bg
        i['hse_metallic'] = hse.metallic
        hse = hse.bg
    if hse is None:
      i['hse_bandgap'] = None
      i['hse_metallic'] = None
    i['hse_location'] = hse_f.replace(hse_fname,'.save')
    hse_lst.append(i)
  return hse_lst


def process_save(fname,counts):
  '''This gets the bandgap, dopant, etc for a range of save directories'''
  f = open(fname).read().split('\n')
  from math import gcd
  string = ''
  cmmon = []
  for i in counts:
    cmmon.append(counts[i])
  if len(cmmon) == 1:
    div = cmmon[0]
  else: 
    div = cmmon[0]
    for c in cmmon[1::]:
        div = gcd(div , c)

  counts_sort = sorted([(i,int(counts[i]/div)) for i in counts],key=lambda x: x[1])
  b_ = [i[0] for i in counts_sort]
  phase = ''
  for i in counts_sort:
    if i[1] == 1:
      phase += i[0]
    else:
      phase += i[0] + str(i[1])

  d = []
  for i in f:
    if len(i) == 0: continue
    try:
      y = QE.Struct()
      y.XML_Process(i)
      x = bg(i.replace('.save','.out'))
      t = [zz for zz in y.atoms if ''.join([tt for tt in zz if not tt.isdigit()]).strip() not in b_]
      p = [zz for zz in y.atoms if ''.join([tt for tt in zz if not tt.isdigit()]).strip() in b_]
      coun = {}
      for zz in p:
        c = ''.join([tt for tt in zz.split() if not tt.isdigit()]) 
        if c in coun: coun[c] += 1
        else: coun[c] = 1
      loc = phase
      loc = b_[0] if coun[b_[0]] < counts[b_[0]] else b_[1]
      if len(t) > 0:
        dopant = ''.join([zz for zz in t[0] if not zz.isdigit()]).strip()
      else:
        dopant = phase
      d.append({'location':i.strip(),'dopant':dopant,'bandgap':x.bg,'metallic':x.metallic,'dop_sub':loc,'phase':phase,'natom':len(y.atoms)})
    except:
      continue
  return phase,d

if __name__ == '__main__':
  phase,d = process_save('save.out',{'Ti':8,'O':16})

  hse_location = '/data/llentz/Charge-Density/HSE/HSE/TiO2/Big/c-len'
  t = get_hse(d,phase,hse_location)
  json.dump(t, open(phase + '.hse.json','w'))

