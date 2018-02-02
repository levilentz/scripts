#!/usr/bin/env python
import numpy as np
import sys

vbm = None
cbm = None
bg = None
class bg:
  def __init__( self, direct):
    self.vbm = None
    self.cbm = None
    self.bg = None
    self.Fermi = None
    self.occu = None
    self.energy = None
    self.metallic = False
    self.get_bandgap(direct)
  
  def __call__ ( self,direct):
    self.get_bandgap(direct)
  
  def get_bandgap(self,direct):
    self.occu = []
    self.energy = []
    f = open(direct)

    for i in f:
      if "End of self-consistent" in i:
        self.occu = []
        self.energy = []
      if ' k =' in i:
        data = []
        i = next(f)
        i = next(f)
        while len(i.split()) > 0:
          for jj in i.split():
            try:
              data.append(float(jj))
            except ValueError:
              jjt = [xx for xx in jj.split('-') if xx]
              for xx in jjt:
                data.append(-float(xx))
          i = next(f)
        self.energy.append(data)
      if 'occupation numbers' in i:
        data = []
        while len(i.split()) > 0:
          i = next(f)
          for jj in i.split():
            data.append(float(jj.strip()))
        self.occu.append(data)
      if 'Fermi' in i:
        self.Fermi = float(i.split()[4])
    f.close() 
    minmax = []
    metallic = False
    if len(self.occu) > 0:
      for c,i in enumerate(self.occu):
        t = np.array(i)
        if 0.0 in t:
          idx = np.where(t==0.0)[0][0]
          te = np.array(self.energy[c])
          minmax.append([te[idx-1],te[idx]])
        else:
          metallic = True 
      if not metallic:
        minmax = np.array(minmax)
        self.vbm = np.max(minmax[:,0])
        self.cbm = np.min(minmax[:,1])
        self.bg = self.cbm - self.vbm
      else:
        self.vbm = None
        self.cbm = None
        self.metallic = True
        self.bg = 0.0
    else:
      minmax = []
      for c,i in enumerate(self.energy):
        t = np.array(i)
        t = np.subtract(i,self.Fermi)
        m_ = np.argwhere(t>0)
        minmax.append([t[m_[0][0]-1],t[m_[0][0]]]) 
      minmax = np.array(minmax)
      self.vbm = np.max(minmax[:,0])
      self.cbm = np.min(minmax[:,1])
      self.bg = self.cbm - self.vbm

if __name__ == "__main__":
  import sys
  x = bg(sys.argv[1])
  print(x.bg)
