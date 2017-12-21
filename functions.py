
import numpy as np
import json
import copy

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

def prophet_list(pname):
  try:
    d_ = open(pname)
    t = []
    for i in d_:
      if '-----------------------' in i:
        i = next(d_)
        while len(i.split()) > 0:
          if len(i.split()) == 5:
            t_ = i.split()
            t.append({'prediction':float(t_[1]),'target':float(t_[2]),'natom':int(t_[3]),'train':t_[4]})
          i = next(d_)
    d_.close()
    return t
  except:
    print('error opening ',pname)
    return 0
       

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

def construct_json(*args):
  t = []
  for i in args:
    t_ = json.load(open(i))
    for j in t_: t.append(j)
  return t

def rot_dir(theta,dir_ = 'x'):
  t = np.pi*theta/180
  s = np.sin
  c = np.cos
  if dir_=='x':
    return np.array([[1,0,0],[0,c(t),-s(t)],[0,s(t),c(t)]])
  elif dir_=='y':
    return np.array([[c(t),0,s(t)],[0,1,0],[-s(t),0,c(t)]])
  elif dir_ == 'z':
    return np.array([[c(t),-s(t),0],[s(t),c(t),0],[0,0,1]])

def rotate_QE(x,theta,dir_='x'):
  y = copy.deepcopy(x)
  rot_matrix = rot_dir(theta,dir_=dir_)
  l = []
  for i in y.lattice:
    l.append(y.lattice[i])
  l = np.transpose(np.array(l))
  l = np.dot(rot_matrix,l)
  for j in y.atoms:
    y.atoms[j] = np.transpose(np.dot(rot_matrix,np.transpose(y.atoms[j])))
  t = ['a','b','c']
  for c,i in enumerate(np.transpose(l)):
    y.lattice[t[c]] = i
  return y

hammett = {'F':	0.34,'NH2':-0.16,'H':0,'COCl':0.51,'CF3':0.43,'OH':0.12,'NHNO2':0.91}
