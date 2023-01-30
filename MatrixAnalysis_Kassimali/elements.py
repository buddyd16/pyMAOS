# -*- coding: utf-8 -*-
'''
BSD 3-Clause License
Copyright (c) 2023, Donald N. Bockoven III
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
import numpy as np

class R2Truss():
    
    def __init__(self, inode, jnode, material, section, uid):
        
        self.type = "TRUSS"
        
        self.inode = inode
        self.jnode = jnode
        
        self.material = material
        self.section = section
        
        self.uid = uid
        
        self.end_forces_local = {}
        
        self.end_forces_global = {}
    
    @property
    def length(self):
        '''
        Sets the member length from the i and j nodes
        '''

        self._length = self.inode.distance(self.jnode)
        
        return self._length
    
    def k(self):
        
        E = self.material.E
        A = self.section.Area
        L = self.length
        
        k = np.matrix([[A*E/L,0,0,-1*A*E/L,0,0],
                   [0,0,0,0,0,0],
                   [0,0,0,0,0,0],
                   [-A*E/L,0,0,A*E/L,0,0],
                   [0,0,0,0,0,0],
                   [0,0,0,0,0,0]])
        
        return k
    
    def T(self):
        
        c = (self.jnode.x-self.inode.x)/self.length
        s = (self.jnode.y-self.inode.y)/self.length
        
        T = np.matrix([[c,s,0,0,0,0],
                       [-s,c,0,0,0,0],
                       [0,0,1,0,0,0],
                       [0,0,0,c,s,0],
                       [0,0,0,-s,c,0],
                       [0,0,0,0,0,1]])
        
        return T
    
    def kglobal(self):
        
        k = self.k()
        T = self.T()
        
        kglobal = np.matmul(np.matmul(np.transpose(T),k),T)
        
        return kglobal
    
    def Flocal(self, loadcase):
        # Gloabl nodal displacement vector
        D = np.zeros(6)
        
        iD = self.inode.displacements[loadcase]
        jD = self.jnode.displacements[loadcase]
        
        # Populate Displacement Vector
        D[0] = iD[0]
        D[1] = iD[1]
        D[2] = iD[2]
        D[3] = jD[0]
        D[4] = jD[1]
        D[5] = jD[2]
        
        dlocal = np.matmul(self.T(),D)

        FL = np.matmul(self.k(),dlocal.T)
        
        self.end_forces_local[loadcase] = FL
        
    def Fglobal(self, loadcase):
        
        # Gloabl nodal displacement vector
        D = np.zeros(6)
        
        iD = self.inode.displacements[loadcase]
        jD = self.jnode.displacements[loadcase]
        
        # Populate Displacement Vector
        D[0] = iD[0]
        D[1] = iD[1]
        D[2] = iD[2]
        D[3] = jD[0]
        D[4] = jD[1]
        D[5] = jD[2]
        
        # global stiffness matrix
        KG = self.kglobal()
        
        FG = np.matmul(KG, D)
        
        self.end_forces_global[loadcase] = FG
        
        self.Flocal(loadcase)

        return FG


class R2Frame():
    
    def __init__(self, inode, jnode, material, section, uid):
        
        self.type = "FRAME"
        
        self.inode = inode
        self.jnode = jnode
        
        self.material = material
        self.section = section
        
        self.uid = uid
    
    @property
    def length(self):
        '''
        Sets the member length from the i and j nodes
        '''

        self._length = self.inode.distance(self.jnode)
        
        return self._length
    
    def k(self):
        
        E = self.material.E
        Ixx = self.section.Ixx
        A = self.section.Area
        L = self.length
        
        k = np.matrix([[A*E/L,0,0,-1*A*E/L,0,0],
                       [0,(12*E*Ixx)/(L*L*L),(6*E*Ixx)/(L*L),0,(-12*E*Ixx)/(L*L*L),(6*E*Ixx)/(L*L)],
                       [0,(6*E*Ixx)/(L*L),(4*E*Ixx)/L,0,(-6*E*Ixx)/(L*L),(2*E*Ixx)/L],
                       [-A*E/L,0,0,A*E/L,0,0],
                       [0,(-12*E*Ixx)/(L*L*L),(-6*E*Ixx)/(L*L),0,(12*E*Ixx)/(L*L*L),(-6*E*Ixx)/(L*L)],
                       [0,(6*E*Ixx)/(L*L),(2*E*Ixx)/L,0,(-6*E*Ixx)/(L*L),(4*E*Ixx)/L]])
        
        return k
    
    def T(self):
        
        c = (self.jnode.x-self.inode.x)/self.length
        s = (self.jnode.y-self.inode.y)/self.length
        
        T = np.matrix([[c,s,0,0,0,0],
                       [-s,c,0,0,0,0],
                       [0,0,1,0,0,0],
                       [0,0,0,c,s,0],
                       [0,0,0,-s,c,0],
                       [0,0,0,0,0,1]])
        
        return T
    
    def kglobal(self):
        
        k = self.k()
        T = self.T()
        
        kglobal = np.matmul(np.matmul(np.transpose(T),k),T)
        
        return kglobal