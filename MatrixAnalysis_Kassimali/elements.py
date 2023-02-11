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
import numpy as np

import loading as loadtypes


class R2Truss:
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
        """
        Sets the member length from the i and j nodes
        """

        self._length = self.inode.distance(self.jnode)

        return self._length

    def k(self):

        E = self.material.E
        A = self.section.Area
        L = self.length

        k = np.matrix(
            [
                [A * E / L, 0, 0, -1 * A * E / L, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [-A * E / L, 0, 0, A * E / L, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
            ]
        )

        return k

    def T(self):

        c = (self.jnode.x - self.inode.x) / self.length
        s = (self.jnode.y - self.inode.y) / self.length

        T = np.matrix(
            [
                [c, s, 0, 0, 0, 0],
                [-s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, s, 0],
                [0, 0, 0, -s, c, 0],
                [0, 0, 0, 0, 0, 1],
            ]
        )

        return T

    def kglobal(self):

        k = self.k()
        T = self.T()

        kglobal = np.matmul(np.matmul(np.transpose(T), k), T)

        return kglobal

    def Dglobal(self, loadcase):
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

        return D

    def Dlocal(self, loadcase):

        Dglobal = self.Dglobal(loadcase)

        Dlocal = np.matmul(self.T(), Dglobal)

        return Dlocal

    def Flocal(self, loadcase):

        Dlocal = self.Dlocal(loadcase)

        FL = np.matmul(self.k(), Dlocal.T)

        self.end_forces_local[loadcase] = FL

    def Fglobal(self, loadcase):

        Dglobal = self.Dglobal(loadcase)

        # global stiffness matrix
        KG = self.kglobal()

        FG = np.matmul(KG, Dglobal)

        self.end_forces_global[loadcase] = FG

        self.Flocal(loadcase)

        return FG


class R2Frame:
    def __init__(self, inode, jnode, material, section, uid):

        self.type = "FRAME"

        self.inode = inode
        self.jnode = jnode

        self.material = material
        self.section = section

        self.uid = uid

        self.loads = []

        self.end_forces_local = {}

        self.end_forces_global = {}

        self.fixed_end_forces = {}

    @property
    def length(self):
        """
        Sets the member length from the i and j nodes
        """

        self._length = self.inode.distance(self.jnode)

        return self._length

    def add_point_load(self, p, a, case, direction, location_percent=False):
        """

        Parameters
        ----------
        p : TYPE
            DESCRIPTION.
        a : TYPE
            DESCRIPTION.
        case : TYPE
            DESCRIPTION.
        direction : TYPE
            DESCRIPTION.
        location_percent : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """

        if location_percent:
            a = (a / 100) * self.length
        if direction == "Y" or direction == "X":
            # Load is applied in the gloabl axis

            c = (self.jnode.x - self.inode.x) / self.length
            s = (self.jnode.y - self.inode.y) / self.length

            if direction == "Y":

                pyy = c * p
                pxx = s * p
            else:
                pyy = -1 * s * p
                pxx = c * p
            self.loads.append(
                loadtypes.R2_Axial_Load(pxx, a, self, loadcase=case)
            )
            self.loads.append(
                loadtypes.R2_Point_Load(pyy, a, self, loadcase=case)
            )
        else:
            # Load is applied in the local member axis

            if direction == "xx":

                self.loads.append(
                    loadtypes.R2_Axial_Load(p, a, self, loadcase=case)
                )
            else:

                self.loads.append(
                    loadtypes.R2_Point_Load(p, a, self, loadcase=case)
                )

    def add_distributed_load(
        self, wi, wj, a, b, case, direction, location_percent=False
    ):
        """

        Parameters
        ----------
        wi : TYPE
            DESCRIPTION.
        wj : TYPE
            DESCRIPTION.
        a : TYPE
            DESCRIPTION.
        b : TYPE
            DESCRIPTION.
        case : TYPE
            DESCRIPTION.
        direction : TYPE
            DESCRIPTION.
        location_percent : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """

        if location_percent:
            a = (a / 100) * self.length
            b = (b / 100) * self.length
        if direction == "Y" or direction == "X":
            # Load is applied in the gloabl axis

            c = (self.jnode.x - self.inode.x) / self.length
            s = (self.jnode.y - self.inode.y) / self.length

            if direction == "Y":

                wyyi = c * wi
                wyyj = c * wj
                wxxi = s * wi
                wxxj = s * wj
            else:
                wyyi = -1 * s * wi
                wyyj = -1 * s * wj
                wxxi = c * wi
                wxxj = c * wj
            self.loads.append(
                loadtypes.R2_Axial_Linear_Load(
                    wxxi, wxxj, a, b, self, loadcase=case
                )
            )
            self.loads.append(
                loadtypes.R2_Linear_Load(wyyi, wyyj, a, b, self, loadcase=case)
            )
        else:
            # Load is applied in the local member axis

            if direction == "xx":

                self.loads.append(
                    loadtypes.R2_Axial_Linear_Load(
                        wi, wj, a, b, self, loadcase=case
                    )
                )
            else:

                self.loads.append(
                    loadtypes.R2_Linear_Load(wi, wj, a, b, self, loadcase=case)
                )

    def add_moment_load(self, m, a, case, location_percent=False):
        """

        Parameters
        ----------
        m : TYPE
            DESCRIPTION.
        a : TYPE
            DESCRIPTION.
        case : TYPE
            DESCRIPTION.
        location_percent : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """

        if location_percent:
            a = (a / 100) * self.length
        self.loads.append(loadtypes.R2_Point_Moment(m, a, self, loadcase=case))

    def FEF(self, case):
        """

        Parameters
        ----------
        case : TYPE
            DESCRIPTION.

        Returns
        -------
        fef : TYPE
            DESCRIPTION.

        """

        fef = np.array([0, 0, 0, 0, 0, 0])

        for load in self.loads:

            if load.loadcase == case:
                loadfef = np.array(load.FEF())
                fef = fef + loadfef
        return fef

    def FEFglobal(self, case):
        """

        Parameters
        ----------
        case : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        fef = np.transpose(self.FEF(case))
        T = self.T()

        return np.matmul(np.transpose(T), fef)

    def k(self):
        """

        Returns
        -------
        k : TYPE
            DESCRIPTION.


        """

        E = self.material.E
        Ixx = self.section.Ixx
        A = self.section.Area
        L = self.length

        k = np.matrix(
            [
                [A * E / L, 0, 0, -1 * A * E / L, 0, 0],
                [
                    0,
                    (12 * E * Ixx) / (L * L * L),
                    (6 * E * Ixx) / (L * L),
                    0,
                    (-12 * E * Ixx) / (L * L * L),
                    (6 * E * Ixx) / (L * L),
                ],
                [
                    0,
                    (6 * E * Ixx) / (L * L),
                    (4 * E * Ixx) / L,
                    0,
                    (-6 * E * Ixx) / (L * L),
                    (2 * E * Ixx) / L,
                ],
                [-A * E / L, 0, 0, A * E / L, 0, 0],
                [
                    0,
                    (-12 * E * Ixx) / (L * L * L),
                    (-6 * E * Ixx) / (L * L),
                    0,
                    (12 * E * Ixx) / (L * L * L),
                    (-6 * E * Ixx) / (L * L),
                ],
                [
                    0,
                    (6 * E * Ixx) / (L * L),
                    (2 * E * Ixx) / L,
                    0,
                    (-6 * E * Ixx) / (L * L),
                    (4 * E * Ixx) / L,
                ],
            ]
        )

        return k

    def T(self):
        """

        Returns
        -------
        T : TYPE
            DESCRIPTION.

        """

        c = (self.jnode.x - self.inode.x) / self.length
        s = (self.jnode.y - self.inode.y) / self.length

        T = np.matrix(
            [
                [c, s, 0, 0, 0, 0],
                [-s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, s, 0],
                [0, 0, 0, -s, c, 0],
                [0, 0, 0, 0, 0, 1],
            ]
        )

        return T

    def kglobal(self):
        """

        Returns
        -------
        kglobal : TYPE
            DESCRIPTION.

        """

        k = self.k()
        T = self.T()

        kglobal = np.matmul(np.matmul(np.transpose(T), k), T)

        return kglobal

    def Dglobal(self, loadcase):
        """

        Parameters
        ----------
        loadcase : TYPE
            DESCRIPTION.

        Returns
        -------
        D : TYPE
            DESCRIPTION.

        """

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

        return D

    def Dlocal(self, loadcase):
        """

        Parameters
        ----------
        loadcase : TYPE
            DESCRIPTION.

        Returns
        -------
        Dlocal : TYPE
            DESCRIPTION.

        """

        Dglobal = self.Dglobal(loadcase)

        Dlocal = np.matmul(self.T(), Dglobal)

        return Dlocal

    def Flocal(self, loadcase):
        """

        Parameters
        ----------
        loadcase : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        Dlocal = self.Dlocal(loadcase)
        Qf = np.reshape(self.FEF(loadcase), (-1, 1))

        FL = np.matmul(self.k(), Dlocal.T)

        print(Qf)
        print(FL)
        print(FL + Qf)

        self.end_forces_local[loadcase] = FL + Qf

    def Fglobal(self, loadcase):
        """

        Parameters
        ----------
        loadcase : TYPE
            DESCRIPTION.

        Returns
        -------
        FG : TYPE
            DESCRIPTION.

        """

        Dglobal = self.Dglobal(loadcase)
        Qfg = self.FEFglobal(loadcase)

        # global stiffness matrix
        KG = self.kglobal()

        FG = np.matmul(KG, Dglobal)

        self.end_forces_global[loadcase] = FG + Qfg

        self.Flocal(loadcase)

        return FG + Qfg
