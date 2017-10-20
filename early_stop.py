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

def to_pd(f,d_,train,df=None):
  if df is not None:
    df_ = df.copy()
    for i in d_:
      df_.loc[i,'target'] = d_[i]['target']
      df_.loc[i,'prediction'] = d_[i]['prophet']
      df_.loc[i,'train'] = train
    return df_
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

def process(fname,df=None,executable='PROPhet',np=32,db=None):
  np = str(np)
  nsave,checkpoint,nint = get_restart(fname)
  valf = convert(fname,'val.dat','FILE')
  out = open('earlystop.out','w')
  len_ = len(open('val.dat').read().split('\n')[:-1])
  np = str(len_) if len_ < 32 else str(32)
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

  t = pm('train.dat.out','train.dat')
  v = pm('val.dat.out','val.dat')
  if df is not None:
    df = to_pd(f,t,'train',df=df)
    df = to_pd(f,v,'val',df=df)
    t = get_net(fname='bfgs_file')
    print(t)
    if db is not None:
      F_pkl = pkl.load(open(db,'rb'))
      F_pkl[os.getcwd()] = {'description':t,'df':df}
    pkl.dump(F_pkl,open(db,'wb'))
    df.to_csv('data.csv')
  else:
    f = open('data.csv','w')
    f.write('target,prediction,train\n')
    to_pd(f,t,'train')
    to_pd(f,v,'train')
    f.close()

def construct_df(j):
  if j is not None:
    _ = pd.read_json(j)
    _.set_index('location',inplace=True)
    return _
  else:
    raise ValueError('json file does not exist',j)
  
if __name__ == "__main__":
  df = construct_df('/data/llentz/codeplayground/data/all.json')
  d = process('bfgs_file',df=df,executable='PROPhet',db='/data/llentz/codeplayground/data/Database.pkl')
