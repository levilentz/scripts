
import numpy as np
import json

def prophet_map(pname,tname):
  '''This is a conversion routine to convert PROPhet out put into a dictionary with the PK being the directory'''
  try:
    d_ = open(tname).read().split('\n')[:-1]
    d = []
    for i in d_: 
      if i is not '':
        d.append(i.split()[0])
    d_ = d
  except:
    print('error opening ',tname)
    return 0
  p_ = {}
  with open(pname,'r') as f:
    for i in f:
      while 'System' not in i: 
        i = next(f)
      if 'train' in i.lower(): train = True
      else: train = False
      i = next(f)
      i = next(f)
      cnt = 0
      while len(i.split()) > 0:
        s_ = i.split()
        if train:
          t_ = {'prophet':s_[1],'target':s_[2],'train':s_[3]}
        else:
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
color_keys = {'H': '#FFFFFF', 'He': '#D9FFFF', 'Li': '#CC80FF', 'Be': '#C2FF00', 'B': '#FFB5B5', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D', 'O2': '#FFAE00', 'F': '#90E050', 'Ne': '#B3E3F5', 'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000', 'S': '#FFFF30', 'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00', 'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'Ti1': '#BFC2C7', 'Ti2': '#BFC2C7', 'V': '#A6A6AB', 'V1': '#A6A6AB', 'V2': '#A6A6AB', 'Cr': '#8A99C7', 'Cr1': '#8A99C7', 'Cr2': '#8A99C7', 'Mn': '#9C7AC7', 'Mn1': '#9C7AC7', 'Mn2': '#9C7AC7', 'Fe': '#FFA800', 'Fe1': '#FFA200', 'Fe2': '#FFD200', 'Co': '#F090A0', 'Co1': '#05004C', 'Co2': '#388786', 'Co3': '#67CAC9', 'Ni': '#50D050', 'Ni1': '#50D050', 'Ni2': '#50D050', 'Cu': '#808080', 'Cu1': '#808080', 'Cu2': '#606060', 'Zn': '#7D80B0', 'Ga': '#C28F8F', 'Ge': '#668F8F', 'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929', 'Kr': '#5CB8D1', 'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0', 'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C', 'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F', 'In': '#A67573', 'Sn': '#668080', 'Sb': '#9E63B5', 'Te': '#D47A00', 'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F', 'Ba': '#00C900', 'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7', 'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7', 'Tb': '#30FFC7', 'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675', 'Tm': '#00D452', 'Yb': '#00BF38', 'Lu': '#00AB24', 'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB', 'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0', 'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00', 'At': '#754F45', 'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00', 'Ac': '#70ABFA', 'Th': '#00BAFF', 'Pa': '#00A1FF', 'U': '#008FFF', 'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2', 'Cm': '#785CE3', 'Bk': '#8A4FE3', 'Cf': '#A136D4', 'Es': '#B31FD4', 'Fm': '#B31FBA', 'Md': '#B30DA6', 'No': '#BD0D87', 'Lr': '#C70066', 'Rf': '#CC0059', 'Db': '#D1004F', 'Sg': '#D90045', 'Bh': '#E00038', 'Hs': '#E6002E', 'Mt': '#EB0026'}
