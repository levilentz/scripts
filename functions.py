
import numpy as np
import json

def prophet_map(pname,tname):
  '''This is a conversion routine to convert PROPhet out put into a dictionary with the PK being the directory'''
  try:
    d_ = open(tname).read().split('\n')[:-1]
    d = []
    for i in d_: 
      if i is not '':
        d.append(i)
    d_ = d
  except:
    print('error opening ',tname)
    return 0
  p_ = {}
  with open(pname,'r') as f:
    for i in f:
      while 'System         Prediction       Target' not in i: 
        i = next(f)
      i = next(f)
      i = next(f)
      cnt = 0
      while len(i.split()) > 0:
        s_ = i.split()
        t_ = {'prophet':s_[1],'target':s_[2]}
        p_[d_[cnt]] = t_
        cnt += 1
        i = next(f)
      break
  return p_

def temperature(scf_file,wd):
  '''old routine for AIMD'''
  f = open(scf_file)
  temp = []
  for i in f:
    if 'temperature' in i and len(i.split()) == 4: 
      temp.append(i.split()[2])
  f.close()
  for c,i in enumerate(temp):
    f_ = open(wd + '/' + str(c+1) + '.save/temperature','w')
    f_.write(i)
    f_.close()

def get_network_info(fname='bfgs_file'):
  '''extracts relevant network information for PROPhet'''
  f = open(fname).read().split('\n')
  d = {}
  for i in f:
    if 'hidden' in i:
      l = i.find('=')
      d['network'] = i[l+1:len(i)].strip()
    if 'downsample' in i:
      l = i.find('=')
      d['downsample'] = i[l+1:len(i)].strip()
    if 'precondition' in i:
      if '1' in i and '#' not in i: 
        d['precondition'] = True
      elif '#' in i:
        d['precondition'] = False
      else:
        d['precondition'] = False
  return d 