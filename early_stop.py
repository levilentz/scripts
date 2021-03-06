#!/usr/bin/env python

import os
from functions import prophet_map as pm
import re
import sqlite3
import pandas as pd
import numpy as np
import os
import pickle as pkl
from functions import get_network_info as get_net

def to_pd(d_,train,df=None,f=None):
  if df is not None:
    print(d_)
    for i in d_:
      df.loc[i,'target'] = d_[i]['target']
      df.loc[i,'prediction'] = d_[i]['prophet']
      if 'train' in list(d_[i].keys()):
        df.loc[i,'train'] = d_[i]['train']
      else:
        df.loc[i,'train'] = train
  else:
    for i in d_:
      f.write(','.join([str(zz) for zz in [d_[i]['target'],d_[i]['prophet'],train]]) + '\n')

def analysis(s):
  d = []
  s_ = s.read().split('\n')
  try:
    for c,i in enumerate(s_):
      if 'System         Prediction       Target' in i:
        c += 2
        i = s_[c]
        while len(i.split()) > 0:
          d.append([float(i.split()[1]),float(i.split()[2])])
          c += 1
          i = s_[c]
        break
    d = np.array(d)
    del_ = np.max(np.abs(d[:,0] - d[:,1]))
    rmse = np.sqrt(np.sum((d[:,0]-d[:,1])**2)/len(d))
  except:
    print(s_)
    raise "error with PROPhet"
  return del_,rmse

def convert(s,include,chkpoint):
  f = open(s)
  inc = False
  chk = False
  val_f = ""
  for i in f:
    if 'checkpoint_in' in i.lower():
      val_f += 'checkpoint_in = ' + chkpoint + '\n'
      chk = True
    elif 'include' in i.lower():
      val_f += 'include = ' + include + '\n'
      inc = True
    else:
      val_f += i
  if not inc:
    val_f += 'include = ' + include + '\n'
  if not chk:
    val_f += 'checkpoint_in = ' + chkpoint + '\n'
  f.close()
  return val_f

def get_restart(s):
  f = open(s)
  nsave = None 
  checkpoint = None
  nint = None 
  for i in f:
    if i[0] == '#': continue
    if 'nsave' in i.lower():
      nsave = i[i.find('=') + 1:].strip()
    if 'checkpoint_out' in i.lower():
      checkpoint = i[i.find('=') + 1:].strip()
    if 'niterations' in i.lower():
      nint = i[i.find('=') + 1:].strip()
  f.close()
  return nsave,checkpoint,nint

def process(fname,bout=None,df=None,executable='PROPhet',np=32,db=None,d=None):
  if bout is not None:
    if d is None: 
      t = []
      f = open(bout)
      for i in f:
        if 'Iteration   ' in i:
          i = next(f)
          i = next(f)
          while len(i.split()) == 4:
            try:
              t.append(i.split())
              i = next(f)
            except: break
      f.close()
      t = sorted(t,key= lambda x: float(x[2]))
      d = [(int(t[0][0]),)]
    nsave,checkpoint,nint = get_restart(fname)
    print(nsave,checkpoint,nint)
    nsave = int(nsave)
    c = int(d[0][0])
    correct = round((c/nsave)+1)*nsave
    valf = convert(fname,'train.dat','FILE')
    if not os.path.isfile(checkpoint + '_' + str(correct)): 
      if os.path.isfile(checkpoint + '_' +  str(int(correct) - int(nsave))):
        correct = correct - nsave
        chkpoint = checkpoint + '_' + str(correct)
      elif int(correct) - int(nsave) == int(nint):
        chkpoint = checkpoint
      else:  
        raise ValueError(correct)
    else: chkpoint = checkpoint + '_' + str(correct)
    f = open('val_temp','w')
    f.write(valf.replace('FILE',chkpoint))
    f.close()
    #t = os.popen('mpirun -np {np} {prop} -in val_temp -validate | tee train.dat.out'.format(prop=executable,np=np))
    t = os.popen('mpirun -np 32 PROPhet -in val_temp -validate | tee train.dat.out'.format(prop=executable,np=np)).read()
    #print(t)
    #funct = open(checkpoint + '_' + str(correct)).read()
    funct = open(chkpoint).read()
    t_file = ['train.dat']
    to_pkl(db=db,df=df,t_file=t_file,funct=funct)
    return
  np = str(np)
  nsave,checkpoint,nint = get_restart(fname)
  valf = convert(fname,'val.dat','FILE')
  len_ = len(open('val.dat').read().split('\n')[:-1])
  np = str(len_) if len_ < 32 else str(32)
  if d is None:
    out = open('earlystop.out','w')
    d = []
    out.write('step,rmse,max\n')
    for i in range(100,int(nint),int(nsave)):
      if not os.path.isfile(checkpoint + '_' + str(i)):
        break
      f = open('val_temp','w')
      f.write(valf.replace('FILE',checkpoint + '_' + str(i)))
      f.close()
      t = os.popen('mpirun -np {np} {prop} -in val_temp -validate'.format(prop=executable,np=np))
      del_,rmse = analysis(t)
      d.append((i,rmse,del_))
      print(i,rmse,del_)
    d = sorted(d,key=lambda x: x[1])
    for i in d:
      out.write(','.join([str(zz) for zz in i]) + '\n')
    out.close()
  f = open('val_temp','w')
  f.write(valf.replace('FILE',checkpoint + '_' + str(d[0][0])))
  f.close()
  t = os.popen('mpirun -np {np} {prop} -in val_temp -validate > val.dat.out'.format(prop=executable,np=np)).read()
  f = open('train_temp','w')
  f.write(valf.replace('FILE',checkpoint + '_' + str(d[0][0])).replace('val.dat','train.dat'))
  f.close()
  t = os.popen('mpirun -np {np} {prop} -in train_temp -validate > train.dat.out'.format(prop=executable,np=np)).read()
  f = open('test_temp','w')
  f.write(valf.replace('FILE',checkpoint + '_' + str(d[0][0])).replace('val.dat','test.dat'))
  f.close()
  t = os.popen('mpirun -np {np} {prop} -in test_temp -validate > test.dat.out'.format(prop=executable,np=np)).read()
  t_file = ['train.dat','val.dat','test.dat']
  #t_out = ['train.dat.out','val.dat.out','test.dat.out']
  #flag = ['train','val','test']
  funct = open(checkpoint + '_' + str(d[0][0])).read()
  to_pkl(db=db,df=df,t_file=t_file,funct=funct)

def to_pkl(db=None,fname='bfgs_file',df=None,t_file=['train.dat'],f=None,funct=None):
  if df is not None:
    for c,i in enumerate(t_file):
      t = pm(i + '.out',i)
      to_pd(t,i.replace('.dat',''),df=df)
    df = df.dropna()
    t = get_net(fname='bfgs_file')
    if db is not None:
      F_pkl = pkl.load(open(db,'rb'))
      F_pkl[os.getcwd()] = {'description':t,'df':df.T.to_dict(),'functional':funct} #storing the dataframe as dict for version control
      pkl.dump(F_pkl,open(db,'wb'))
    df.to_csv('data.csv')
  else:
    f.write('target,prediction,train\n')
    to_pd(t,'train',f=f)
    to_pd(v,'train',f=f)

def construct_df(j):
  if j is not None:
    _ = pd.read_json(j)
    _.set_index('location',inplace=True)
    _['target'] = None
    _['prediction'] = None
    _['train'] = None
    return _
  else:
    raise ValueError('json file does not exist',j)

def split_val(df,val_file='val.dat'):
  t = open('val.dat').split('\n')[:-1]
  d = df.ix[t]
  val_temp = pd.DataFrame()
  test_temp = pd.DataFrame()
  for i in d.phase.unique():
    for j in d[d.phase == i].dopant.unique():
      v_t = d[(d.phase == i) & (d.dopant == j)]
      vt = v_t.sample(frac=0.5)
      tt = v_t.drop(vt.index)
      val_temp = val_temp.append(vt)
      test_temp = test_temp.append(tt)
  f = open('val.dat','w')
  for i in val_temp.index:
    f.write(i + '\n')
  f.close()
  f = open('test.dat','w')
  for i in test_temp.index:
    f.write(i + '\n')
  f.close()
  
if __name__ == "__main__":
  #df = construct_df('/data/llentz/Charge-Density/no_Phosphate/data/all.json')
  #d = process('bfgs_file',df=df,executable='PROPhet',db='/data/llentz/codeplayground/data/Database.pkl')
  df = construct_df('/data/llentz/Charge-Density/HSE/data/all.hse.json')
  d = process('bfgs_file',df=df,executable='PROPhet',db='/data/llentz/Charge-Density/HSE/data/database.hse.pkl',bout='train.bfgs')
