
import numpy as np
import pandas as pd
import json
import copy

def get_pdos(key,directory = '.'):
  d = []
  cnt = []
  files = glob.glob(directory + '/' + '*(' + key + ')*')
  for i in files:
    lb = i.find('#') + 1
    rb = i.find('(',lb)
    if i[lb:rb] not in cnt: cnt.append(i[lb:rb])
  for i in [files[0]]:
    f = open(i)
    next(f)
    next(f)
    d_ = []
    for jj in f:
      d_.append([float(zz) for zz in [jj.split()[0],jj.split()[1]]])
    f.close()
    if len(d) == 0: d = np.array(d_)
    else:
      d[:,1] += np.array(d_)[:,1]
  return d.tolist(),len(cnt)

def get_files(directory = '.'):
  files = glob.glob(directory + '/' + '*pdos_atm*')
  pdos = {}
  cnt = {}
  for i in files:
    lb = i.find('(')
    rb = i.find(')',lb)
    key = i[lb+1:rb]
    if key not in pdos:
      pdos[key],cnt[key] = get_pdos(key,directory=directory)
  return pdos,cnt

def get_fermi(fname='vc-relax.out',directory = '.'):
  f = open(directory + '/' + fname)
  fermi = []
  for i in f:
    if 'Fermi' in i:
      fermi.append(float(i.split()[-2]))
  return fermi

def prophet_map(pname,tname):
  '''This is a conversion routine to convert PROPhet out put into a dictionary with the PK being the directory'''
  try:
    d_ = open(tname).read().split('\n')[:-1]
    d = []
    train = []
    for i in d_: 
      if i is not '':
        i_ = i.split()
        d.append(i_[0])
        if len(i_) > 1:
          train.append(i_[1])
        else:
          train.append(i_[0])
    d_ = d
  except:
    print('error opening ',tname)
    return 0
  p_ = {}
  with open(pname,'r') as f:
    for i in f:
      while 'System' not in i: 
        i = next(f)
      i = next(f)
      i = next(f)
      cnt = 0
      while len(i.split()) > 0:
        if 'warning' in i.lower(): continue
        s_ = i.split()
        t_ = {'prophet':float(s_[1]),'target':float(s_[2]),'train':train[cnt]}
        p_[d_[cnt]] = t_
        cnt += 1
        i = next(f)
      break
  return p_

def subscript(string):
  t = ''
  for i in string:
    if i.isdigit(): t+= '$_' + i + '$'
    else: t += i
  return t

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

def split_df(df,**kwargs):
  '''This will add a column to a dataframe and split it into train, val, and test datasets'''
  if not kwargs:
    kwargs = **{'frac':0.80}
  df['train'] = None
  train = df.sample(**kwargs)
  df.loc[train.index,'train'] = 'train'
  rem = df.drop(train.index)
  test = rem.sample(frac=0.50)
  df.loc(test.index,'train') = 'test'
  df.loc(rem.drop(test.index).index,'train') = 'val'

def upf(stru):
  atom,cell = stru.return_params()
  t = {'H': 'H 1.0079 H.upf', 'He': 'He 4.0026 He.upf', 'Li': 'Li 6.941 Li.upf', 'Be': 'Be 9.0122 Be.upf', 'B': 'B 10.811 B.upf', 'C': 'C 12.0107 C.upf', 'N': 'N 14.0067 N.upf', 'O': 'O 15.9994 O.upf', 'F': 'F 18.9984 F.upf', 'Ne': 'Ne 20.1797 Ne.upf', 'Na': 'Na 22.9897 Na.upf', 'Mg': 'Mg 24.305 Mg.upf', 'Al': 'Al 26.9815 Al.upf', 'Si': 'Si 28.0855 Si.upf', 'P': 'P 30.9738 P.upf', 'S': 'S 32.065 S.upf', 'Cl': 'Cl 35.453 Cl.upf', 'Ar': 'Ar 39.948 Ar.upf', 'K': 'K 39.0983 K.upf', 'Ca': 'Ca 40.078 Ca.upf', 'Sc': 'Sc 44.9559 Sc.upf', 'Ti': 'Ti 47.867 Ti.upf', 'V': 'V 50.9415 V.upf', 'Cr': 'Cr 51.9961 Cr.upf', 'Mn': 'Mn 54.938 Mn.upf', 'Fe': 'Fe 55.845 Fe.upf', 'Co': 'Co 58.9332 Co.upf', 'Ni': 'Ni 58.6934 Ni.upf', 'Cu': 'Cu 63.54600000000001 Cu.upf', 'Zn': 'Zn 65.39 Zn.upf', 'Ga': 'Ga 69.723 Ga.upf', 'Ge': 'Ge 72.64 Ge.upf', 'As': 'As 74.9216 As.upf', 'Se': 'Se 78.96 Se.upf', 'Br': 'Br 79.904 Br.upf', 'Kr': 'Kr 83.8 Kr.upf', 'Rb': 'Rb 85.4678 Rb.upf', 'Sr': 'Sr 87.62 Sr.upf', 'Y': 'Y 88.9059 Y.upf', 'Zr': 'Zr 91.22399999999999 Zr.upf', 'Nb': 'Nb 92.9064 Nb.upf', 'Mo': 'Mo 95.94 Mo.upf', 'Tc': 'Tc 98.0 Tc.upf', 'Ru': 'Ru 101.07 Ru.upf', 'Rh': 'Rh 102.9055 Rh.upf', 'Pd': 'Pd 106.42 Pd.upf', 'Ag': 'Ag 107.8682 Ag.upf', 'Cd': 'Cd 112.411 Cd.upf', 'In': 'In 114.818 In.upf', 'Sn': 'Sn 118.71 Sn.upf', 'Sb': 'Sb 121.76 Sb.upf', 'Te': 'Te 127.6 Te.upf', 'I': 'I 126.9045 I.upf', 'Xe': 'Xe 131.293 Xe.upf', 'Cs': 'Cs 132.9055 Cs.upf', 'Ba': 'Ba 137.327 Ba.upf', 'La': 'La 138.9055 La.upf', 'Ce': 'Ce 140.116 Ce.upf', 'Pr': 'Pr 140.9077 Pr.upf', 'Nd': 'Nd 144.24 Nd.upf', 'Pm': 'Pm 145.0 Pm.upf', 'Sm': 'Sm 150.36 Sm.upf', 'Eu': 'Eu 151.964 Eu.upf', 'Gd': 'Gd 157.25 Gd.upf', 'Tb': 'Tb 158.9253 Tb.upf', 'Dy': 'Dy 162.5 Dy.upf', 'Ho': 'Ho 164.9303 Ho.upf', 'Er': 'Er 167.25900000000001 Er.upf', 'Tm': 'Tm 168.9342 Tm.upf', 'Yb': 'Yb 173.04 Yb.upf', 'Lu': 'Lu 174.967 Lu.upf', 'Hf': 'Hf 178.49 Hf.upf', 'Ta': 'Ta 180.9479 Ta.upf', 'W': 'W 183.84 W.upf', 'Re': 'Re 186.207 Re.upf', 'Os': 'Os 190.23 Os.upf', 'Ir': 'Ir 192.217 Ir.upf', 'Pt': 'Pt 195.078 Pt.upf', 'Au': 'Au 196.9665 Au.upf', 'Hg': 'Hg 200.59 Hg.upf', 'Tl': 'Tl 204.3833 Tl.upf', 'Pb': 'Pb 207.2 Pb.upf', 'Bi': 'Bi 208.9804 Bi.upf', 'Po': 'Po 209.0 Po.upf', 'At': 'At 210.0 At.upf', 'Rn': 'Rn 222.0 Rn.upf', 'Fr': 'Fr 223.0 Fr.upf', 'Ra': 'Ra 226.0 Ra.upf', 'Ac': 'Ac 227.0 Ac.upf', 'Th': 'Th 232.0381 Th.upf', 'Pa': 'Pa 231.0359 Pa.upf', 'U': 'U 238.0289 U.upf', 'Np': 'Np 237.0 Np.upf', 'Pu': 'Pu 244.0 Pu.upf', 'Am': 'Am 243.0 Am.upf'}
  atm = {}
  for i in atom.split('\n')[:-1]:
    _ = i.split()[0]
    if _ in atm: continue
    atm[_] = t[_]
  _ = '\n'.join(atm[i] for i in atm)
  return _


hammett = {'F': 0.34,'NH2':-0.16,'H':0,'COCl':0.51,'CF3':0.43,'OH':0.12,'NHNO2':0.91}
color_keys = {'H': '#FFFFFF', 'He': '#D9FFFF', 'Li': '#CC80FF', 'Be': '#C2FF00', 'B': '#FFB5B5', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D', 'O2': '#FFAE00', 'F': '#90E050', 'Ne': '#B3E3F5', 'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000', 'S': '#FFFF30', 'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00', 'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'Ti1': '#BFC2C7', 'Ti2': '#BFC2C7', 'V': '#A6A6AB', 'V1': '#A6A6AB', 'V2': '#A6A6AB', 'Cr': '#8A99C7', 'Cr1': '#8A99C7', 'Cr2': '#8A99C7', 'Mn': '#9C7AC7', 'Mn1': '#9C7AC7', 'Mn2': '#9C7AC7', 'Fe': '#FFA800', 'Fe1': '#FFA200', 'Fe2': '#FFD200', 'Co': '#F090A0', 'Co1': '#05004C', 'Co2': '#388786', 'Co3': '#67CAC9', 'Ni': '#50D050', 'Ni1': '#50D050', 'Ni2': '#50D050', 'Cu': '#808080', 'Cu1': '#808080', 'Cu2': '#606060', 'Zn': '#7D80B0', 'Ga': '#C28F8F', 'Ge': '#668F8F', 'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929', 'Kr': '#5CB8D1', 'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0', 'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C', 'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F', 'In': '#A67573', 'Sn': '#668080', 'Sb': '#9E63B5', 'Te': '#D47A00', 'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F', 'Ba': '#00C900', 'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7', 'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7', 'Tb': '#30FFC7', 'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675', 'Tm': '#00D452', 'Yb': '#00BF38', 'Lu': '#00AB24', 'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB', 'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0', 'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00', 'At': '#754F45', 'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00', 'Ac': '#70ABFA', 'Th': '#00BAFF', 'Pa': '#00A1FF', 'U': '#008FFF', 'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2', 'Cm': '#785CE3', 'Bk': '#8A4FE3', 'Cf': '#A136D4', 'Es': '#B31FD4', 'Fm': '#B31FBA', 'Md': '#B30DA6', 'No': '#BD0D87', 'Lr': '#C70066', 'Rf': '#CC0059', 'Db': '#D1004F', 'Sg': '#D90045', 'Bh': '#E00038', 'Hs': '#E6002E', 'Mt': '#EB0026'}
