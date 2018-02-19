#!/usr/bin/env python
import sys, os
import numpy as np
from numpy import cos as cos
from numpy import sin as sin
import json
from collections import Counter
import sqlite3 as lite
from xml.etree import ElementTree
from lxml import etree as ET
import copy

class index():
  def __init__(self):
    self.keys = {}

  def sanitize(self,s):
    return ''.join([i for i in s if not i.isdigit()]).strip()

  def key(self,specie):
    if specie not in self.keys:
      self.keys[specie] = 0
    else:
      self.keys[specie] += 1
    return specie + str(self.keys[specie])

  def reset(self):
    self.keys = {}

class Struct:
  def __init__(self,direct = None):
    self.BOHRtoA = 0.529177249
    self.RYtoeV = 13.605698066
    self.program = "QE"
    self.version = ""
    self.volume = 0.0
    self.alat = 0.0
    self.natoms = 0
    self.nat = 0
    self.nelect = 0
    self.Ecut = 0.0
    self.RhoCut = 0.0
    self.Econv = 0.0
    self.Exch = ""
    self.energy = 0.0
    self.natoms = 0
    self.bandgap = 0.0
    self.bands = 0
    self.lattice = {'a':np.zeros(3),'b':np.zeros(3),'c':np.zeros(3)}
    self.atoms = {}
    self.norms = {'a':0.0,'b':0.0,'c':0.0}
    self.angles = {'alpha':0.0,'beta':0.0,'gamma':0.0}
    self.kpts = 0
    self.bnddiagram = False
    self.FermiTest = False
    self.Fermi = 0.0
    self.JSON = ""
    self.email = ""
    self.atomindex = index()
    self.noband = False
    if direct is not None:
      if os.path.isfile(direct):
        self.File_Process(direct)
      elif os.path.isdir(direct):
        self.XML_Process(direct)

  def From_Crystal(self):
    #vol = np.sqrt(1 - np.cos(self.angles['alpha']*np.pi/180.)**2 - np.cos(self.angles['beta']*np.pi/180.)**2 - np.cos(self.angles['gamma']*np.pi/180.)**2 + 2*np.cos(self.angles['alpha']*np.pi/180.)*np.cos(self.angles['beta']*np.pi/180.)*np.cos(self.angles['gamma']*np.pi/180.))
    #temp1 = [self.norms['a'], self.norms['b']*np.cos(self.angles['gamma']*np.pi/180.), self.norms['c']*np.cos(self.angles['beta']*np.pi/180.)]
    #temp2 = [0,self.norms['b']*np.sin(self.angles['gamma']*np.pi/180.),self.norms['c']*((np.cos(self.angles['alpha']*np.pi/180.)-np.cos(self.angles['beta']*np.pi/180.)*np.cos(self.angles['beta']*np.pi/180.))/np.sin(self.angles['gamma']*np.pi/180.))]
    #temp3 = [0,0,self.norms['c']*vol/np.sin(self.angles['gamma']*np.pi/180.)]
    #conversion = np.vstack((np.asarray(temp1),np.asarray(temp2),np.asarray(temp3))) 
    t = []
    for i in ['a','b','c']:
      t.append(list(self.lattice[i]))
    return np.transpose(np.array(t))
    #return np.linalg.inv(np.transpose(np.array(t)))
  
  def to_Crystal(self):
    self.Normalize()
    conversion = np.linalg.inv(self.From_Crystal())
    print(conversion)
    print("ATOM_POSITIONS {crystal}")
    for i in self.atoms:
      #tmp = copy.deepcopy(self.atoms[i])
      tmp = np.dot(conversion,self.atoms[i])
      print(self.atomindex.sanitize(i) + " " + ' '.join([str(round(j,9)) for j in tmp])) 

  def wrap_Cell(self):
    self.Normalize()
    conversion = np.linalg.inv(self.From_Crystal())
    c_1 = self.From_Crystal()
    for i in self.atoms:
      tmp = np.dot(conversion,self.atoms[i])
      for c,j in enumerate(tmp):
        if tmp[c] < 0.0: 
          while tmp[c] <= 0.0:
            tmp[c] += 1
        elif tmp[c] > 1.0: 
          while tmp[c] >= 1.0:
            tmp[c] -= 1
      self.atoms[i] = np.dot(c_1,tmp)
      
    

  def RDF(self,rcut=5.0,dr=0.1):
    supcell = []
    rho = self.natoms/self.volume
    radius = np.arange(0,rcut + dr,dr)
    R = {}
    for i in radius:
      R[i] = 0.0
    conversion = np.linalg.inv(self.From_Crystal())
    crystal = self.From_Crystal()
    max_trans = {}
    for i in self.norms:
      max_trans[i] = int(np.ceil(1.1*rcut/self.norms[i]))
    trans_index = []
    for i in range(-max_trans['a'],max_trans['a']+1):
      for j in range(-max_trans['b'],max_trans['b']+1):
        for k in range(-max_trans['c'],max_trans['c']+1):
          trans_index.append([i,j,k])
    y = copy.deepcopy(self)
    for i in y.atoms:
      y.atoms[i] = np.dot(conversion,y.atoms[i])
    
    for i in y.atoms:
      pos1 = copy.copy(y.atoms[i])
      for j in y.atoms:
        for z in trans_index:
          pos2 = np.add(y.atoms[j],z)
          delta = np.linalg.norm(np.dot(crystal,np.subtract(pos2,pos1)))
          if delta <= rcut:
            for ZZ in range(0,len(radius)-1):
              if radius[ZZ] < delta < radius[ZZ+1]:
                R[radius[ZZ]] += 1.0/(rho*4.*np.pi*(radius[ZZ]**2)*dr)
    data = []
    for i in R:
      data.append([i,R[i]/(float(self.natoms))])
    return data
    
    
  def CIF(self,filename=""):
    self.Normalize()
    ciffile = "data_global\n"
    ciffile += "_chemical_name " + self.to_Formula().strip() + "\n"
    for i in ['a','b','c']:
      ciffile += '_cell_length_' + i.strip() + " " + str(self.norms[i]) + "\n"
    for i in ['alpha','beta','gamma']:
      ciffile += '_cell_angle_' + i.strip() + " " + str(self.angles[i]) + "\n"
    ciffile += '_cell_volume_ ' + str(self.volume) + "\n"
    ciffile += "loop_\n"
    ciffile += "_atom_site_label\n"
    ciffile += "_atom_site_fract_x\n"
    ciffile += "_atom_site_fract_y\n"
    ciffile += "_atom_site_fract_z\n"
    conversion = np.linalg.inv(self.From_Crystal())
    counter = 0
    for i in self.atoms:
      tmp = copy.deepcopy(self.atoms[i])
      dot = np.dot(conversion,tmp)
      ciffile += self.atomindex.sanitize(i) + " " + " ".join([str(x) for x in dot]) + "\n"
    if filename == "":
      print(ciffile)
    else:
      f = open(filename, 'w')
      f.write(ciffile)
      f.close()
  def to_XYZ(self,filename = None):
    printstring = ''
    printstring += str(len(self.atoms)) + '\n\n'
    for i in self.atoms:
      tmp = [round(x,5) for x in self.atoms[i]]
      string = self.atomindex.sanitize(i) + " " + " ".join([str(x) for x in self.atoms[i]])
      printstring += string + "\n"
    if filename == None:
      print(printstring)
    else:
      f = open(filename,'w')
      f.write(printstring)
      f.close()

  def return_params(self):
    atoms = ''
    for i in self.atoms:
      tmp = [round(x,5) for x in self.atoms[i]]
      string = self.atomindex.sanitize(i) + " " + " ".join([str(x) for x in self.atoms[i]])
      atoms += string + "\n"
    cell = ''
    for i in ['a','b','c']:
      string = " ".join([str(round(x,5)) for x in self.lattice[i]])
      cell += string + "\n"
    return atoms,cell

  def print(self):
    printstring = "ATOMIC_POSITIONS {angstrom}\n"
    for i in self.atoms:
      tmp = [round(x,5) for x in self.atoms[i]]
      string = self.atomindex.sanitize(i) + " " + " ".join([str(x) for x in self.atoms[i]])
      printstring += string + "\n"
    printstring += "CELL_PARAMETERS {angstrom}\n"
    for i in ['a','b','c']:
      string = " ".join([str(round(x,5)) for x in self.lattice[i]])
      printstring += string + "\n"
    return printstring

  def XML_Process(self,dirstring): #Need to talk to Alexie About this, does not store total energy
    try:
      f = open(dirstring + "/data-file.xml")
    except:
      raise ValueError("Cannot open " +  dirstring + '/data-file.xml')
    tree = ElementTree.parse(f)
    f.close()
    MAP = {'a':'a1','b':'a2','c':'a3'}
    try:
      self.version = tree.find('./HEADER/CREATOR').attrib['VERSION']
      self.nat = int(tree.find('./IONS/NUMBER_OF_SPECIES').text)
      self.atomindex.reset()
      self.Exch = tree.find('./EXCHANGE_CORRELATION/DFT').text.rstrip().lstrip()
      self.Nelec = float(tree.find('./BAND_STRUCTURE_INFO/NUMBER_OF_ELECTRONS').text.strip())
      self.kpts = float(tree.find('./BAND_STRUCTURE_INFO/NUMBER_OF_K-POINTS').text)
      self.fermi = 27.2114*float(tree.find('./BAND_STRUCTURE_INFO/FERMI_ENERGY').text)
      self.alat = self.BOHRtoA*float(tree.find('./CELL/LATTICE_PARAMETER').text)
      self.Ecut = 27.2114*float(tree.find('./PLANE_WAVES/WFC_CUTOFF').text)
      self.RhoCut = 27.2114*float(tree.find('./PLANE_WAVES/WFC_CUTOFF').text)
      #self.beta = 0.0
    except:
      self.version = None
      self.natoms = None
      self.nat = None
      self.atomindex.reset()
      self.Exch = None
      self.Nelec = None
      self.kpts = None
      self.fermi = None
      self.alat = None
      self.Ecut = None
      self.RhoCut = None
    try:
      self.energy = float(tree.find('./TOTAL_ENERGY').text)
    except:
      self.energy = None
    for i in ['a','b','c']:
      try:
        for c,j in enumerate(tree.find('./CELL/DIRECT_LATTICE_VECTORS/' + MAP[i]).text.split()):
          self.lattice[i][c] = self.BOHRtoA*float(j)
      except:
        for c,j in enumerate(tree.find('./DIRECT_LATTICE_VECTORS/' + MAP[i]).text.split()):
          self.lattice[i][c] = self.BOHRtoA*float(j)
    self.natoms = int(tree.find('./IONS/NUMBER_OF_ATOMS').text)
    for i in range(0,self.natoms):
      tmp = tree.find('./IONS/ATOM.' + str(i + 1)).attrib
      species = tmp['SPECIES']
      array = [ self.BOHRtoA*float(j) for j in tmp['tau'].split()]
      self.atoms[self.atomindex.key(species)] = np.array([array[0],array[1],array[2]])
    self.Normalize()
    #test = tree.find('./IONS/ATOM.1')
    #print(test)

  def to_Latex(self):
    self.Normalize()
    t = '''\\begin{table}[]
\centering
\caption{CAPTION}
\label{my-label}
\\bigskip
\\begin{tabular}{lr}
\hline
 Parameter & Value  \\\\
\hline
 PARAMS
\hline
\end{tabular}
\end{table}'''
    params = ''
    for c,i in enumerate(self.norms):
      params += i.upper() + ' & ' + str(round(self.norms[i],2)) + ' \\\\' + '\n'
    for i in self.angles:
      params += '$\\' + i + '$ & ' + str(round(self.angles[i],2)) + ' \\\\' + '\n'
    print(t.replace('PARAMS',params).replace('CAPTION',self.to_Formula() + ' unit cell parameters. Length units are in \\AA and angles are in degrees.'))
    t1 = '''\\begin{table}[]
\centering
\caption{CAPTION}
\\bigskip
\label{my-label}
\\begin{tabular}{lrrr}
\hline
 Atom & X & Y & Z  \\\\
\hline
 PARAMS
\hline
\end{tabular}
\end{table}'''
    params_1 = ''
    conversion = np.linalg.inv(self.From_Crystal())
    print()
    for i in self.atoms:
      #tmp = copy.deepcopy(self.atoms[i])
      tmp = np.dot(conversion,self.atoms[i])
      params_1 += self.atomindex.sanitize(i) + " & " + ' & '.join([str(round(j,3)) for j in tmp]) + '\\\\ \n'
    print(t1.replace('PARAMS',params_1).replace('CAPTION',self.to_Formula() + ' atomic positions. Here positions are shown in fractional coordinates'))

  def File_Process(self,filestring):
    try:
      f = open(filestring,'r')
    except:
      print("Cannot open %s" % filestring)
      #sys.exit(1)
      raise IOError
    linenum = 0
    for i in f:
      i = i.lower()
      #if "ERROR" in i.upper():
        #print("There is an error in this calculation")
        #sys.exit(2)
      if linenum < 1000:
        if 'nat' in i.lower():
          self.natoms = int(''.join(zz for zz in i.strip() if zz.isdigit()))
        if "lattice parameter (a_0)" in i:
          self.alat = float(i.split()[5])
        if "number of k points=" in i:
          self.kpts = int(i.split()[4])
          next
        if "Program PWSCF" in i:
          self.version = i.split()[2].replace('v.','')
          next
        if "lattice parameter (alat)" in i:
          self.alat = float(i.split()[4])*self.BOHRtoA
          next
        if "number of Kohn-Sham states" in i:
          self.bands = int(i.split()[4]) 
        if "unit-cell volume" in i and "new" not in i:
          self.volume = float(i.split()[3])*(self.BOHRtoA**3.0)
          next
        if "number of atoms/cell" in i:
          self.natoms = int(i.split()[4])
          next
        if "number of atomic types" in i:
          self.nat = int(i.split()[5])
          next
        if "number of electrons" in i:
          self.nelect = float(i.split()[4])
          next
        if "kinetic-energy cutoff" in i:
          self.Ecut = float(i.split()[3])*self.RYtoeV
          next
        if "charge density cutoff" in i:
          self.RhoCut = float(i.split()[4])*self.RYtoeV
          next
        if "convergence threshold" in i:
          if len(i.split()) < 4:
            self.Econv = float(i.split()[3])
            next
        if "Exchange-correlation" in i:
          self.Exch = i[i.find('=') + 1:i.find('(')].rstrip()
          next
        if "a(1) =" in i:
          tmp = i.replace('a(1)','').replace('(','').replace('=','').replace(',','').replace(')','').split()
          for j in range(0,3):
            self.lattice['a'][j] = self.alat*float(tmp[j])
          next
        if "a(2) =" in i:
          tmp = i.replace('a(2)','').replace('(','').replace('=','').replace(',','').replace(')','').split()
          for j in range(0,3):
            self.lattice['b'][j] = self.alat*float(tmp[j])
          next
        if "a(3) =" in i:
          tmp = i.replace('a(3)','').replace('(','').replace('=','').replace(',','').replace(')','').split()
          for j in range(0,3):
            self.lattice['c'][j] = self.alat*float(tmp[j])
          next
        if "site n.     atom                  positions (alat units)" in i:
          self.atomindex.reset()
          for j in range(0,self.natoms):
            line = next(f).split()
            self.atoms[self.atomindex.key(line[1])] = np.multiply(np.array([float(line[6]),float(line[7]),float(line[8])]),self.alat)
          next
        if 'nat' in i.lower():
          self.natoms = int(''.join([zz for zz in i if zz.isdigit()]))
          next
      if "!" in i and "ENERGY" in i.upper():
        self.energy= float(i.split()[4])*self.RYtoeV
      if "new unit-cell volume" in i:
        self.volume = float(i.split()[4])*(self.BOHRtoA**3)
      if "cell_parameters" in i:
        if "angstrom" in i:
          for j in ['a','b','c']:
            line = next(f)
            tmp = line.split()
            for k in range(0,3):
              self.lattice[j][k] = float(tmp[k])
        else:
          for j in ['a','b','c']:
            line = next(f)
            tmp = line.split()
            for k in range(0,3):
              self.lattice[j][k] = self.alat*float(tmp[k])
        self.Normalize()
      if "ATOMIC_POSITIONS" in i.upper():
        self.atomindex.reset()
        if "angstrom" in i:
          for j in range(0,self.natoms):
            line = next(f).split()
            self.atoms[self.atomindex.key(line[0])] = np.array([float(line[1]),float(line[2]),float(line[3])])
        if "alat" in i:
          for j in range(0,self.natoms):
            line = next(f).split()
            self.atoms[self.atomindex.key(line[0])] = np.array([self.alat*float(line[1]),self.alat*float(line[2]),self.alat*float(line[3])])
        if "crystal" in i.lower():
          conversion = self.From_Crystal()
          for j in range(0,self.natoms):
            line = next(f).split()
            tmp = np.transpose(np.array([float(line[1]),float(line[2]),float(line[3])]))
            ncoords = np.dot(conversion,tmp)
            self.atoms[self.atomindex.key(line[0])] = np.array([float(ncoords[0]),float(ncoords[1]),float(ncoords[2])])
      if "End of self-consistent calculation" in i:
        if np.floor(self.bands/8.)*8. <= self.bands:
          numlines = int(np.floor(self.bands/8.) + 1)
          remainder = int(self.bands - np.floor(self.bands/8.)*8.)
        else: 
          numlines = int(np.floor(self.bands/8.))
          remainder = 0
        self.bnddiagram = np.zeros((self.kpts,self.bands))
        counter = 0
        self.noband = False
        while counter < self.kpts:
          line = next(f)
          if "Number of k-points >=" in line: 
            self.noband = True
            break
          if "k =" in line:
            line = next(f)
            counter1 = 0
            for j in range(0,numlines):
              line = next(f)
              '''
              for k in range(0,len(line.split())):
                self.bnddiagram[counter][counter1 + k] = float(line.split()[k])
              '''
              counter1 += 8
            counter += 1  
        next
      if "highest occupied, lowest unoccupied level (ev)" in i:
        self.bandgap = float(i.split()[7]) - float(i.split()[6])
        next
      if "the Fermi energy is" in i:
        self.Fermi = float(i.split()[4])
        self.FermiTest = True
        next
      linenum += 1
    f.close()
    self.Normalize()
    self.noband = True
    if self.FermiTest == True and self.noband == False:
      self.bnddiagram = np.subtract(self.bnddiagram,self.Fermi)
      emin = np.zeros(self.kpts)
      emax = np.zeros(self.kpts)
      counter = 0
      for j in self.bnddiagram:
        emin[counter] = j[np.searchsorted(j,  0.0,side='right')-1]
        emax[counter] = j[np.searchsorted(j,  0.0,side='right')]
        counter += 1
      self.bandgap = float(np.min(emax-emin))

  def to_JSON(self):
    if self.JSON == "":
      for i in self.lattice:
        self.lattice[i] = self.lattice[i].tolist() 
      self.bnddiagram = self.bnddiagram.tolist()
      self.JSON = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=1)

  def to_Formula(self):
    from math import gcd
    string = ''
    cmmon = []
    for i in self.atomindex.keys:
      cmmon.append(self.atomindex.keys[i] + 1)
    if len(cmmon) == 1:
      div = cmmon[0]
    else: 
      div = cmmon[0]
      for c in cmmon[1::]:
          div = gcd(div , c)
    for i in self.atomindex.keys:
      t_ = str(int((self.atomindex.keys[i]+1)/div))
      if t_ == '1': t_ = ''
      string += i + t_
    return string

  def to_database(self):
    con = None
    try:
      con = lite.connect('QE.db')
      cur = con.cursor()    
      insert_command = 'INSERT INTO QE(A,B,C,ALPHA,BETA,GAMMA,VOLUME,NATOMS,FORMULA,BANDGAP,ENERGY,PROGRAM,VERSION) VALUES ('
      for i in('a','b','c'):
        insert_command += str(self.norms[i]) + ','
      for i in('alpha','beta','gamma'):
        insert_command += str(self.angles[i]) + ','
      insert_command += str(self.volume) + ',' +  str(self.natoms) + ',"' + self.to_Formula() + '",' + str(self.bandgap) + ',' + str(self.energy) + ',"QE",' + '"' + self.version + '");'
      cur.execute(insert_command)
      lid = cur.lastrowid
      JSON = self.to_JSON()
      insert_command = "INSERT INTO RAW_DATA VALUES (" + str(lid) + ",'" +  self.JSON + "');"
      cur.execute(insert_command)
      self.to_File(lid)
    except lite.Error:
      print("Unable to insert into the database")
      #sys.exit(3)
      raise IOError
    finally:
      if con:
        con.commit()
        con.close()

  def Normalize(self):
    try:
      for i in ['a','b','c']:
        self.norms[i] = np.linalg.norm(self.lattice[i])
      self.angles['alpha'] = np.arccos(np.dot(self.lattice['b'],self.lattice['c'])/(self.norms['c']*self.norms['b'])) * 180./np.pi
      self.angles['gamma'] = np.arccos(np.dot(self.lattice['a'],self.lattice['b'])/(self.norms['a']*self.norms['b'])) * 180./np.pi
      self.angles['beta'] = np.arccos(np.dot(self.lattice['a'],self.lattice['c'])/(self.norms['a']*self.norms['c'])) * 180./np.pi
      self.volume = np.dot(self.lattice['a'],np.cross(self.lattice['b'],self.lattice['c']))
    except:
      print("Lattice undefined")
      raise ValueError
      #sys.exit(4)

  def to_File(self,lid):
    with open(str(lid) + '.json','w') as f:
      f.write(self.JSON)
  
  def to_XML(self,fname):
    root = ET.Element("Root") 
    DIRECT = ET.SubElement(root,'DIRECT_LATTICE_VECTORS')
    UNITS = ET.SubElement(DIRECT,'UNITS_FOR_DIRECT_LATTICE_VECTORS')
    UNITS.set('UNITS',"Bohr")
    lattice = {}
    convert = {'a1':'a','a2':'b','a3':'c'}
    for i in range(1,4):
      key = 'a' + str(i)
      lattice[key] = ET.SubElement(DIRECT,key)
      lattice[key].set('type','real')
      lattice[key].set('size','3')
      lattice[key].set('columns','3')
      text = ' '.join([str(1.88973*x) for x in self.lattice[convert[key]]]) 
      lattice[key].text = text
    IONS = ET.SubElement(root,'IONS')
    NA = ET.SubElement(IONS,'NUMBER_OF_ATOMS')
    NA.set('type','integer')
    NA.set('size','1')
    NA.text=str(self.natoms)
    NA = ET.SubElement(IONS,'NUMBER_OF_SPECIES')
    NA.set('type','integer')
    NA.set('size','1')
    NA.text = str(self.nat)
    NA = ET.SubElement(IONS,'UNITS_FOR_ATOMIC_POSITIONS')
    NA.set('UNITS','bohr')
    index = 1
    tun = {}
    counter = 1
    for i in self.atoms:
      if self.atomindex.sanitize(i) not in tun:
        tun[self.atomindex.sanitize(i)] = str(index)
        index += 1
      ATOM = ET.SubElement(IONS,'ATOM.' + str(counter))
      ATOM.set('SPECIES',self.atomindex.sanitize(i) + " ")
      ATOM.set('INDEX',tun[self.atomindex.sanitize(i)])
      text = ' '.join([str(x*1.88973) for x in self.atoms[i]])
      ATOM.set('tau',text)
      ATOM.set('if_pos',"1 1 1")
      counter += 1
    TE = ET.SubElement(root,'TOTAL_ENERGY')
    TE.set('UNITS','eV')
    TE.text = str(self.energy)
    f = ET.ElementTree(root)
    f.write(fname + '.xml',pretty_print=True)

  def to_Supercell(self,array,symm=False):
    if isinstance(array,list):
      tmp = copy.deepcopy(self)
      conversion = np.linalg.inv(tmp.From_Crystal())
      for i in tmp.atoms:
        dot = np.dot(conversion,tmp.atoms[i])
        tmp.atoms[i] = dot
      COORDs = copy.copy(tmp.atoms)
      if symm:
        r1 = range(-array[0],array[0] + 1)
        r2 = range(-array[1],array[1] + 1)
        r3 = range(-array[2],array[2] + 1)
      else:
        r1 = range(0,array[0])
        r2 = range(0,array[1])
        r3 = range(0,array[2])
      for i in tmp.atoms:
        for j in r1:
          for k in r2:
            for z in r3:
              if (j == 0) and (k == 0) and (z == 0): #We already have the (0,0,0) structure
                next
              else:
                COORDs[tmp.atomindex.key(tmp.atomindex.sanitize(i))] = np.add(tmp.atoms[i],np.array([j,k,z]))
      for i in COORDs:
        COORDs[i] = np.dot(tmp.From_Crystal(),COORDs[i])
      tmpmap = {'a':0,'b':1,'c':2}
      array = [float(xx) for xx in array]
      for i in tmp.lattice:
        #tmp.lattice[i] = np.array([tmp.lattice[i][0]*float(array[tmpmap[i]]), tmp.lattice[i][1]*float(array[tmpmap[i]]), tmp.lattice[i][2]*float(array[tmpmap[i]])])
        tmp.lattice[i] = tmp.lattice[i]*array[tmpmap[i]]
      tmp.atoms = COORDs
      tmp.natoms = tmp.natoms*(array[0]*array[1]*array[2])
      indexing = ['a','b','c']
      for i in indexing:
        tmp.norms[i] *= array[indexing.index(i)]
      tmp.volume = tmp.volume*(array[0]*array[1]*array[2]) 
      tmp.energy = tmp.energy*(array[0]*array[1]*array[2])
      return tmp
    else:
      print("Invalid supercell dimensions")
  
  def __str__(self):
    return self.print()
  
  def __eq__(self,other):
    if type(other) == type(self):
      diff = []
      keys_ = []
      l = sorted(self.atoms.items(),key=lambda x: (x[1][0],x[1][1],x[1][2]))
      r = sorted(other.atoms.items(),key=lambda x: (x[1][0],x[1][1],x[1][2]))
      for c,i in enumerate(l):
        diff.append(np.linalg.norm(i[1] - r[c][1]))
        keys_.append(self.atomindex.sanitize(i[0]) == self.atomindex.sanitize(r[c][0]))
      m = np.max(np.abs(diff))
      if m > 1e-3: return False
      if False in keys_: return False
      diff = []
      for c,i in enumerate(self.lattice): 
        diff.append(np.linalg.norm(self.lattice[i] - other.lattice[i]))
      m = np.max(np.abs(diff))
      if m > 1e-3: return False
    else: 
      return False
    return True

def main(command):
  test = Struct()
  test.email = command[1]
  if os.path.isfile(command[0]):
    test.File_Process(command[0])
  elif os.path.isdir(command[0]):
    print("Is Dir//process xml here")
    test.XML_Process(command[0])
    print(test.lattice)
  if "@" not in test.email:
    print("Invalid Email Supplied")
    #sys.exit(5)
    raise ValueError
  test.to_database()
  test = None

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Incorrect number of arguments, run as ./QE.py QEOUTPUT_FILE EMAIL")
    sys.exit(6)
  command = [sys.argv[1],sys.argv[2]]
  
  main(command)  
