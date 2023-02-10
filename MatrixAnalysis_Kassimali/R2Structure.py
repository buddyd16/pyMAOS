# -*- coding: utf-8 -*-
"""
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
"""

# from nodes import R2Node
# from elements import R2Truss, R2Frame
# from material import LinearElasticMaterial as Material
# from section import Section
# import loading as loads

import numpy as np


class R2Structure:
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
        self.NDOF = (self.NJD * self.NJ) - self.NR

        # Freedom Map
        self.FM = self.freedom_map()

        # Structure Stiffness Matrix
        self.KSTRUCT = self.Kstructure()

        # Data Stores
        self._D = {}  # Structure Displacement Vector Dictionary

    def freedom_map(self):
        # Freedom Map
        FM = np.zeros(self.NJD * self.NJ)

        # Loop through the nodes mapping free and restrained joint displacements to
        # the Freedom Map (FM). This will facilitate generating the global stiffness
        # matrix in partitioned form.

        j = 0  # starting index for the first free displacement
        k = self.NDOF  # starting index for the first restraint

        for i, node in enumerate(self.nodes):

            for r, restraint in enumerate(node.restraints):

                fmindex = (i) * self.NJD + r

                if restraint == 0:

                    FM[fmindex] = j
                    j += 1
                else:
                    FM[fmindex] = k
                    k += 1
        return FM

    def Kstructure(self):
        """
        Build the structure stiffness matrix orgnized into paritioned form 
        using the freedom map to reposition nodal DOFs

        Returns
        -------
        KSTRUCT: Numpy Matrix
            Structure Stiffness Matrix.

        """

        # Structure Stiffness Matrix
        KSTRUCT = np.zeros([self.NJD * self.NJ, self.NJD * self.NJ])

        for member in self.members:

            # Member global stiffness matrix
            kmglobal = member.kglobal()

            # Freedom map for i and j nodes
            imap = [
                int(self.FM[(member.inode.uid - 1) * self.NJD + r])
                for r in range(self.NJD)
            ]
            imap.extend(
                [
                    int(self.FM[(member.jnode.uid - 1) * self.NJD + r])
                    for r in range(self.NJD)
                ]
            )

            for i in range(self.NJD):

                for y in range(self.NJD):

                    KSTRUCT[imap[i], imap[y]] += kmglobal[i, y]
                    KSTRUCT[imap[i + self.NJD], imap[y]] += kmglobal[i + self.NJD, y]
                    KSTRUCT[imap[i], imap[y + self.NJD]] += kmglobal[i, y + self.NJD]
                    KSTRUCT[imap[i + self.NJD], imap[y + self.NJD]] += kmglobal[
                        i + self.NJD, y + self.NJD
                    ]
        return KSTRUCT

    def nodalforcevector(self, loadcase):
        """
        Build the structure nodal force vector mapped to the same partitions
        as KSTRUCT using the freedom map (FM).

        Returns
        -------
        FG : Numpy Array
            Structure Nodal Force Vector.

        """
        FG = np.zeros(self.NJD * self.NJ)

        for node in self.nodes:

            loads = node.loads.get(loadcase, [0, 0, 0])

            for i, load in enumerate(loads):

                fmindex = (node.uid - 1) * self.NJD + i

                FG[int(self.FM[fmindex])] += load
        return FG

    def solve_linear_static(self, loadcase):
        """
        Perform a linear static solution of the model using the Kff
        and FGf paritions

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        self._verify_stable()

        if self._unstable:

            return 0
        else:

            FG = self.nodalforcevector(loadcase)
            # Slice out the Kff partition from the global structure stiffness
            # Matrix
            self.Kupper = self.KSTRUCT[0 : self.NDOF, 0 : self.NDOF]

            # Slice out the FGf partition from the global nodal force vector
            self.FGupper = FG[0 : self.NDOF]

            # Use Numpy linear Algebra solve function to solve for the
            # displacements at the free nodes.
            U = np.linalg.solve(self.Kupper, self.FGupper)

            # Full Displacement Vector
            # Result is still mapped to DOF via FM
            USTRUCT = np.zeros(self.NJD * self.NJ)

            # Add the resulting free displacements to the appropriate spots in
            # the Full displacement vector
            USTRUCT += np.pad(U, (0, self.NJD * self.NJ - np.shape(U)[0]))

            # store displacement results to the current case to the nodes
            for node in self.nodes:

                uxindex = int(self.FM[(node.uid - 1) * self.NJD + 0])
                uyindex = int(self.FM[(node.uid - 1) * self.NJD + 1])
                rzindex = int(self.FM[(node.uid - 1) * self.NJD + 2])

                node_displacements = [
                    USTRUCT[uxindex],
                    USTRUCT[uyindex],
                    USTRUCT[rzindex],
                ]

                node.displacements[loadcase] = node_displacements
            # compute reactions
            self.compute_reactions(loadcase)

            return U

    def compute_reactions(self, loadcase):
        # Compute Reactions
        for node in self.nodes:

            NL = node.loads.get(loadcase, [0, 0, 0])
            rx = -1 * NL[0]
            ry = -1 * NL[1]
            mz = -1 * NL[2]

            for member in self.members:

                member_FG = member.Fglobal(loadcase)

                if member.inode == node:
                    rx += member_FG[0, 0]
                    ry += member_FG[0, 1]
                    mz += member_FG[0, 2]
                if member.jnode == node:
                    rx += member_FG[0, 3]
                    ry += member_FG[0, 4]
                    mz += member_FG[0, 5]
            node.reactions[loadcase] = [rx, ry, mz]

    def _verify_stable(self):
        """
        Check the diagonal terms of the stiffness matrix against support conditions
        If diagonal term is 0 and the node is unsupported for that DOF then the
        Kmatrix is singular and unstable.

        Returns
        -------
        If unstable returns a dictionary of unstable nodes and degree of freedom marked unstable.

        """

        unstablenodes = []

        for node in self.nodes:

            # Check each DOF of node:
            for i, dof in enumerate(node.restraints):

                fmindex = (node.uid - 1) * self.NJD + i
                val = self.FM[fmindex]

                # value the diagonal position in the stiffness matrix
                kval = self.KSTRUCT[int(val), int(val)]

                if kval == 0 and dof != 1:

                    self._unstable = True
                    unstablenodes.append(
                        f"Node {node.uid} : Unstable for {node.restraints_key[i]}"
                    )
        # add unstable messages to ERROR list
        self._ERRORS.extend(unstablenodes)

        return unstablenodes
