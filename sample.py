#!/usr/bin/env python

import pandas as pd
import json
import numpy as np

def user_input(x):
  with open(x.location + '/user_input','w') as f:
    f.write(str(x.bandgap) + '\n')
    f.write(str(np.log(x.bandgap)) + '\n')
    f.write(str(x.bandgap) + '\n')
d = pd.read_json('all.json')
d = d[(d.metallic == False) & (d.bandgap > 0.15)]
d = d.drop(d[(d.phase == 'GaAs') & (d.dopant != 'GaAs')].index)
d.apply(user_input,axis=1)

#d = pd.read_csv('data.csv')
t_total = pd.DataFrame()
v_total = pd.DataFrame()
test_total = pd.DataFrame()
phase = d['phase'].unique()
for j in phase:
  t_ = d[d.phase == j]
  dopants = t_['dopant'].unique()
  for i in dopants:
    t = t_[t_.dopant == i]
    if len(t) < 10: continue
    train = t.sample(frac=0.8)
    _ = t.drop(train.index)
    val = _.sample(frac=0.50)
    test = _.drop(val.index)
    t_total = t_total.append(train)
    v_total = v_total.append(val)
    test_total = test_total.append(test)

train = open('train.dat','w')
for i in t_total.location.as_matrix():
  train.write(i) 
  train.write('\n')
train.close()
val = open('val.dat','w')
for i in v_total.location.as_matrix():
  val.write(i) 
  val.write('\n')
val.close()
test = open('test.dat','w')
for i in test_total.location.as_matrix():
  test.write(i)
  test.write('\n')
test.close()
