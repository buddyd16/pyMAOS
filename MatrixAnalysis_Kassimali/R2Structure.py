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

from nodes import R2Node
from elements import R2Truss, R2Frame
from material import LinearElasticMaterial as Material
from section import Section
import loading as loads

import matplotlib.pyplot as plt

import numpy as np

class R2Structure():
    
    def __init__(self, nodes, members):
        
        self.nodes = nodes
        self.members = members
        
        # Flags 
        self._unstable = False
        self._Kgenerated = False
        self._ERRORS = []
        
        # Structure Type
        # 2D Structure
        # Ux,Uy, and Rz = 3
        # Number of possible Joint Displacements
        self.NJD = 3
        
        # Number of Joints
        self.NJ = len(self.nodes)
        
        # Number of Members
        self.NM = len(self.members)
        
        # Number of Restraints
        self.NR = sum([sum(n.restraints) for n in self.nodes])
        
        # Degrees of Freedom
        self.NDOF = (self.NJD*self.NJ)-self.NR
        
        # Freedom Map
        self.FM = self.freedom_map()
        
        # Structure Stiffness Matrix
        self.KSTRUCT = self.Kstructure()
        
        # Data Stores
        self._D = {} # Structure Displacement Vector Dictionary


    def freedom_map(self):
        # Freedom Map
        FM = np.zeros(self.NJD*self.NJ)
        
        # Loop through the nodes mapping free and restrained joint displacements to
        # the Freedom Map (FM). This will facilitate generating the global stiffness
        # matrix in partitioned form.
        
        j = 0 # starting index for the first free displacement
        k = self.NDOF # starting index for the first restraint
        
        for i, node in enumerate(self.nodes):
               
            for r, restraint in enumerate(node.restraints):
                
                fmindex = (i)*self.NJD + r
                
                if restraint == 0:
                    
                    FM[fmindex] = j
                    j += 1
                
                else:
                    FM[fmindex] = k
                    k += 1

        return FM


    def Kstructure(self):
        '''
        Build the structure stiffness matrix orgnized into paritioned form 
        using the freedom map to reposition nodal DOFs

        Returns
        -------
        KSTRUCT: Numpy Matrix
            Structure Stiffness Matrix.

        '''
        
        # Structure Stiffness Matrix
        KSTRUCT = np.zeros([self.NJD*self.NJ,self.NJD*self.NJ])
        
        for member in self.members:
            
            # Member global stiffness matrix
            kmglobal = member.kglobal()
        
            # Freedom map for i and j nodes
            imap = [int(self.FM[(member.inode.uid-1)*self.NJD+r]) for r in range(self.NJD)]
            imap.extend([int(self.FM[(member.jnode.uid-1)*self.NJD+r]) for r in range(self.NJD)])
            
            for i in range(self.NJD):
                
                for y in range(self.NJD):
                    
                    KSTRUCT[imap[i],imap[y]] += kmglobal[i, y]
                    KSTRUCT[imap[i+self.NJD],imap[y]] += kmglobal[i+self.NJD,y]
                    KSTRUCT[imap[i],imap[y+self.NJD]] += kmglobal[i,y+self.NJD]
                    KSTRUCT[imap[i+self.NJD],imap[y+self.NJD]] += kmglobal[i+self.NJD,y+self.NJD]
        
        return KSTRUCT
        
        
    def nodalforcevector(self, loadcase):
        '''
        Build the structure nodal force vector mapped to the same partitions
        as KSTRUCT using the freedom map (FM).

        Returns
        -------
        FG : Numpy Array
            Structure Nodal Force Vector.

        '''
        FG = np.zeros(self.NJD*self.NJ)
        
        for node in self.nodes:
            
            loads = node.loads.get(loadcase,[0,0,0])
            
            for i,load in enumerate(loads):
                
                fmindex = (node.uid-1)*self.NJD + i
                
                FG[int(self.FM[fmindex])] += load
        
        return FG


    def solve_linear_static(self, loadcase):
        '''
        Perform a linear static solution of the model using the Kff
        and FGf paritions

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        
        self._verify_stable()
        
        if self._unstable:
            
           return 0
       
        else:
            
            FG = self.nodalforcevector(loadcase)
            # Slice out the Kff partition from the global structure stiffness
            # Matrix
            self.Kupper = self.KSTRUCT[0:self.NDOF,0:self.NDOF]
            
            # Slice out the FGf partition from the global nodal force vector
            self.FGupper = FG[0:self.NDOF]
            
            # Use Numpy linear Algebra solve function to solve for the 
            # displacements at the free nodes.
            U = np.linalg.solve(self.Kupper, self.FGupper)
            
            # Full Displacement Vector
            # Result is still mapped to DOF via FM
            USTRUCT = np.zeros(self.NJD*self.NJ)
            
            # Add the resulting free displacements to the appropriate spots in
            # the Full displacement vector
            USTRUCT += np.pad(U, (0,self.NJD*self.NJ-np.shape(U)[0]))
            
            # store displacement results to the current case to the nodes
            for node in self.nodes:
                
                uxindex = int(self.FM[(node.uid-1)*self.NJD + 0])
                uyindex = int(self.FM[(node.uid-1)*self.NJD + 1])
                rzindex = int(self.FM[(node.uid-1)*self.NJD + 2])
                
                node_displacements = [USTRUCT[uxindex],USTRUCT[uyindex],USTRUCT[rzindex]]
                
                node.displacements[loadcase] = node_displacements
            
            # compute reactions
            self.compute_reactions(loadcase)
            
            return U
    
    def compute_reactions(self, loadcase):
        # Compute Reactions
        for node in self.nodes:
            
            NL = node.loads.get(loadcase,[0,0,0])
            rx = -1*NL[0]
            ry = -1*NL[1]
            mz = -1*NL[2]
            
            for member in self.members:
                
                member_FG = member.Fglobal(loadcase)
                
                if member.inode == node:
                    rx += member_FG[0,0]
                    ry += member_FG[0,1]
                    mz += member_FG[0,2]
                
                if member.jnode == node:
                    rx += member_FG[0,3]
                    ry += member_FG[0,4]
                    mz += member_FG[0,5]
        
            node.reactions[loadcase] = [rx,ry,mz]
        
    def _verify_stable(self):
        '''
        Check the diagonal terms of the stiffness matrix against support conditions
        If diagonal term is 0 and the node is unsupported for that DOF then the
        Kmatrix is singular and unstable.

        Returns
        -------
        If unstable returns a dictionary of unstable nodes and degree of freedom marked unstable.

        '''
        
        unstablenodes = []
        
        for node in self.nodes:
            
            # Check each DOF of node:
            for i, dof in enumerate(node.restraints):
                
                fmindex = (node.uid - 1)*self.NJD + i
                val = self.FM[fmindex]
                
                # value the diagonal position in the stiffness matrix
                kval = self.KSTRUCT[int(val),int(val)]
                
                if kval==0 and dof != 1:
                    
                    self._unstable = True
                    unstablenodes.append(f"Node {node.uid} : Unstable for {node.restraints_key[i]}")
        
        # add unstable messages to ERROR list
        self._ERRORS.extend(unstablenodes)
         
        return unstablenodes


##########################################
##                                      ##
## CONSISTENT UNIT SYTEM REQUIRED !!!!! ##
##                                      ##
##########################################



# Material and Section
M1 = Material(29000)
M2 = Material(10000)
S1 = Section(8,30.8)
S2 = Section(12,0)
S3 = Section(16,0)

# Nodes
N1 = R2Node(0, 0, 1)
N2 = R2Node(288, 0, 2)
N3 = R2Node(576, 0, 3)
N4 = R2Node(864, 0, 4)
N5 = R2Node(288, 216, 5)
N6 = R2Node(576, 216, 6)

# Restraints
N1.restraints = [1,1,1]
N2.restraints = [0,0,1]
N3.restraints = [0,1,1]
N4.restraints = [0,1,1]
N5.restraints = [0,0,1]
N6.restraints = [0,0,1]

# Loading
loadcase = "DL"

N2.loads[loadcase] = [0,-75,0]
N5.loads[loadcase] = [25,0,0]
N6.loads[loadcase] = [0,-60,0]

# Node List
nodes = [N1, N2, N3, N4, N5, N6]


# Elements
T1 = R2Truss(N1, N2, M1, S1, 1)
T2 = R2Truss(N2, N3, M1, S1, 2)
T3 = R2Truss(N3, N4, M2, S3, 3)
T4 = R2Truss(N5, N6, M1, S1, 4)
T5 = R2Truss(N2, N5, M1, S1, 5)
T6 = R2Truss(N3, N6, M1, S1, 6)
T7 = R2Truss(N1, N5, M1, S2, 7)
T8 = R2Truss(N2, N6, M1, S2, 8)
T9 = R2Truss(N3, N5, M1, S2, 9)
T10 = R2Truss(N4, N6, M2, S3, 10)

# T1 = R2Frame(N1, N2, M1, S1, 1)
# T2 = R2Frame(N2, N3, M1, S1, 2)
# T3 = R2Frame(N3, N4, M2, S3, 3)
# T4 = R2Frame(N5, N6, M1, S1, 4)
# T5 = R2Frame(N2, N5, M1, S1, 5)
# T6 = R2Frame(N3, N6, M1, S1, 6)
# T7 = R2Frame(N1, N5, M1, S2, 7)
# T8 = R2Frame(N2, N6, M1, S2, 8)
# T9 = R2Frame(N3, N5, M1, S2, 9)
# T10 = R2Frame(N4, N6, M2, S3, 10)

# Member List
members = [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10]

# Create the 2D Structure
Structure = R2Structure(nodes, members)

Errors = Structure._ERRORS
K = Structure.Kstructure()
U = Structure.solve_linear_static("DL")

print(Errors)
print("Displacements:")
for node in nodes:
    tx = node.displacements["DL"]
    print(f"N{node.uid} -- Ux: {tx[0]:.4E}  Uy:{tx[1]:.4E}  Rz:{tx[2]:.4E}")

print("-"*100)
print("Reactions:")
for node in nodes:
    rx = node.reactions["DL"]
    print(f"N{node.uid} -- Rx: {rx[0]:.4E}  Ry:{rx[1]:.4E}  Mz:{rx[2]:.4E}")

print("-"*100)
print("Member Forces:")
for member in members:
    fx = member.end_forces_local["DL"]
    print(f"M{member.uid} -- Axial: {fx[0,0]:.4E}")

#plot the structure
# fig, ax = plt.subplots()

# for node in nodes:
#     ax.plot(node.x,node.y, marker=".", markersize=20, color='red')
#     ax.plot(node.x_displaced("DL",100),node.y_displaced("DL",100), marker=".", markersize=10, color='gray')

# for member in members:
#     ax.plot([member.inode.x,member.jnode.x],[member.inode.y,member.jnode.y],linewidth=1, color='blue')
#     ax.plot([member.inode.x_displaced("DL",100),member.jnode.x_displaced("DL",100)],[member.inode.y_displaced("DL",100),member.jnode.y_displaced("DL",100)],linewidth=1, color='gray')
    
# ax.grid(True)
# fig.tight_layout()

# plt.show()

# Load test bed
N1 = R2Node(0, 0, 1)
N2 = R2Node(10, 0, 2)

M5 = Material(4176000.0)
S5 = Section(0.02056,0.0014853395061728396)

B1 = R2Frame(N1, N2, M5, S5, 1)

lineload = loads.R2_Axial_Linear_Load(-10, -5, 2, 8, B1)
ptload = loads.R2_Axial_Load(24, 0, B1)
ptload2 = loads.R2_Axial_Load(21, B1.length, B1)

print(lineload.Rix)
print(lineload.Rjx)
print(ptload.Rix)
print(ptload.Rjx)
print(ptload2.Rix)
print(ptload2.Rjx)
print(lineload.Ax)
print(lineload.Dx)
print(ptload.Ax)
print(ptload.Dx)
print(ptload2.Ax)
print(ptload2.Dx)
print(lineload.FEF())
print(ptload.FEF())
print(ptload2.FEF())


A = ptload2.Ax.combine(ptload.Ax.combine(lineload.Ax, 1, 1),1,1)
D = ptload2.Dx.combine(ptload.Dx.combine(lineload.Dx, 12, 12),12,1)

num_stations = 100
eta = [0+i*(1/num_stations)*B1.length for i in range(num_stations+1)]
v = [A.evaluate(i) for i in eta]
fig, ax = plt.subplots()
ax.plot(eta,v, marker=".", markersize=2, color='red')
ax.grid(True)
fig.tight_layout()

plt.show()

