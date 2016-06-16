from __future__ import print_function
import numpy as np

class Cube:
    def __init__(self, filename):
        self.comment = []
        self.sign    = None # sign on number atoms determines if we are dealing
                            # with total density or a molecular orbital.
        self.natoms  = None
        self.origin  = None 
        self.nx      = None 
        self.ny      = None 
        self.nz      = None 
        self.x       = None 
        self.y       = None 
        self.z       = None 
        self.atoms   = []
        with open(filename,'r') as f:
           # Save comments (first two lines)
           for i in xrange(2):
               self.comment.append(f.readline().strip('\n')) 
           # Grab number of atoms and the cubefile origin
           line = f.readline().split()
           self.natoms = abs(int(line[0]))
           self.sign   = np.sign(int(line[0]))
           self.origin = np.array([float(line[1]),
                                   float(line[2]),
                                   float(line[3])])
           # Grab number voxels and axis vector
           for vector in [['nx','x'],['ny','y'],['nz','z']]:
               line = f.readline().split()
               self.__dict__[vector[0]] = int(line[0])
               self.__dict__[vector[1]] = np.array([float(line[1]),
                                                    float(line[2]),
                                                    float(line[3])])
           # Grab atoms and coordinates 
           for atom in xrange(self.natoms):
               line = f.readline().split()
               self.atoms.append([line[0],line[1],line[2],line[3],line[4]])
           # Skip next line ... not sure what it means in new cubegen
           if self.sign < 0:
               f.readline()
           # Grab volumetric data, which finishes out our file
           vals = [float(v) for s in f for v in s.split()]
           # now we need to see what kind of data we have
           # (real/complex?,two-component)
           if len(vals) == self.nx*self.ny*self.nz:
               self.volRA = np.zeros((self.nx,self.ny,self.nz))
               for idx,v in enumerate(vals):
                   self.volRA[idx/(self.ny*self.nz),
                           (idx/self.nz)%self.ny,
                            idx%self.nz] = float(v)
           elif len(vals) == 2*self.nx*self.ny*self.nz:
               # This is the complex one component case
               self.volRA = np.zeros((self.nx,self.ny,self.nz))
               self.volIA = np.zeros((self.nx,self.ny,self.nz))
               for idx,v in enumerate(vals):
                   i = int(np.floor(idx/2))
                   if idx % 2 == 0:
                       self.volRA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 2 == 1:
                       self.volIA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
           elif len(vals) == 4*self.nx*self.ny*self.nz:
               # This is the complex 2-component case, so we have to split the
               # volumetric data into four parts. The density is stored real
               # alpha, imaginary alpha, real beta, imaginary beta 
               self.volRA = np.zeros((self.nx,self.ny,self.nz))
               self.volIA = np.zeros((self.nx,self.ny,self.nz))
               self.volRB = np.zeros((self.nx,self.ny,self.nz))
               self.volIB = np.zeros((self.nx,self.ny,self.nz))
               for idx,v in enumerate(vals):
                   i = int(np.floor(idx/4))
                   if idx % 4 == 0:
                       self.volRA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 4 == 1:
                       self.volIA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 4 == 2:
                       self.volRB[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 4 == 3:
                       self.volIB[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
           else:
               raise NameError, "cube file not valid"

    def write_out(self,filename,data='RA'):
        # Generate new cube file
        # String variable 'data' specifies which component you want
        #  RA = real alpha
        #  IA = imag alpha
        #  RB = real beta 
        #  IB = imag beta
        if data == 'RA':
            volume = self.volRA
        elif data == 'IA':
            volume = self.volIA
        elif data == 'RB':
            volume = self.volRB
        elif data == 'IB':
            volume = self.volIB

        with open(filename, 'w') as f:
            for i in self.comment:
                print(str(i),file=f)
            print(" %4d %.6f %.6f %.6f" % (self.sign*self.natoms, 
                self.origin[0], self.origin[1],self.origin[2]),file=f)
            print(" %4d %.6f %.6f %.6f" % (self.nx, self.x[0], self.x[1],
                self.x[2]), file=f)
            print(" %4d %.6f %.6f %.6f" % (self.ny, self.y[0], self.y[1],
                self.y[2]), file=f)
            print(" %4d %.6f %.6f %.6f" % (self.nz, self.z[0], self.z[1],
                self.z[2]), file=f)
            for atom in self.atoms:
                print(" %s %s %s %s %s" % (atom[0], atom[1], atom[2], atom[3],
                    atom[4]), file=f)
            if self.sign < 0:
                print("    1    "+str(self.natoms),file=f)
            for ix in xrange(self.nx):
                for iy in xrange(self.ny):
                    for iz in xrange(self.nz):
                        print(" %.5e " % volume[ix,iy,iz],end="",file=f)
                        if (iz % 6 == 5):
                            print('',file=f)
                    print('',file=f)
       

if __name__ == '__main__':
    atom = Cube('twoc.cube')
    atom.write_out('twoc_ra.cube',data='RA')
    atom.write_out('twoc_ia.cube',data='IA')
    atom.write_out('twoc_rb.cube',data='RB')
    atom.write_out('twoc_ib.cube',data='IB')
    
