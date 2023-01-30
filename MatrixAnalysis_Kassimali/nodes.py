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

import math

class R2Node():
    
    def __init__(self, x, y, uid):
        '''
        
        Parameters
        ----------
        x : float
            x-axis coordinate of the node in R2.
        y : float
            y-axis coordinate of the node in R2.
        uid : int
            unique node number.

        Returns
        -------
        None.

        '''
        
        self.uid = uid
        self.x = x
        self.y = y
        
        # Restraints [Ux, Uy, Rz]
        self.restraints_key = ["Ux","Uy", "Rz"]
        self.restraints = [0,0,0]
        
        # Dict of Loads by case
        self.loads = {}
        
        # Dict of Displacements by case
        self.displacements = {}
        
        # Dict of Reactions by case
        self.reactions = {}
        
    def __str__(self):
        
        str = f"Node:{self.uid}\n"
        str += f"({self.x},{self.y})\n"
        
        if sum(self.restraints) != 0:
            str += "Restraints\n"
            str += "-"*11+"\n"
            
            for i,r in enumerate(self.restraints):
                if r != 0:
                    str += f"{self.restraints_key[i]}"
        return str
    
    def x_displaced(self, loadcase, scale=1.0):
        
        delta = self.displacements.get(loadcase,[0,0,0])
        
        return self.x + (delta[0]*scale)
    
    
    def y_displaced(self, loadcase, scale=1.0):
        
        delta = self.displacements.get(loadcase,[0,0,0])
        
        return self.y + (delta[1]*scale)
    
    
    def distance(self, other):
        '''
        
        Parameters
        ----------
        other : R2Node
            another node defined by this class.

        Returns
        -------
        distance: float
            Euclidean distance in R2

        '''
        
        dx = (self.x-other.x)
        dy = (self.y-other.y)
        
        return math.sqrt((dx*dx)+(dy*dy))