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

    def Dglobal(self, combo):
        # Gloabl nodal displacement vector
        D = np.zeros(6)

        iD = self.inode.displacements[combo]
        jD = self.jnode.displacements[combo]

        # Populate Displacement Vector
        D[0] = iD[0]
        D[1] = iD[1]
        D[2] = iD[2]
        D[3] = jD[0]
        D[4] = jD[1]
        D[5] = jD[2]

        return D

    def Dlocal(self, combo):

        Dglobal = self.Dglobal(combo)

        Dlocal = np.matmul(self.T(), Dglobal)

        return Dlocal

    def Flocal(self, combo):

        Dlocal = self.Dlocal(combo)

        FL = np.matmul(self.k(), Dlocal.T)

        self.end_forces_local[combo] = FL

    def Fglobal(self, combo):

        Dglobal = self.Dglobal(combo)

        # global stiffness matrix
        KG = self.kglobal()

        FG = np.matmul(KG, Dglobal)

        self.end_forces_global[combo] = FG

        self.Flocal(combo)

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

        # Flags
        self._stations = False
        self._loaded = False

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

        self._stations = False
        self._loaded = True

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

        self._stations = False
        self._loaded = True

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

        self._stations = False
        self._loaded = True

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

    def Dglobal(self, combo):
        """

        Parameters
        ----------
        combo : TYPE
            DESCRIPTION.

        Returns
        -------
        D : TYPE
            DESCRIPTION.

        """

        # Gloabl nodal displacement vector
        D = np.zeros(6)

        iD = self.inode.displacements[combo]
        jD = self.jnode.displacements[combo]

        # Populate Displacement Vector
        D[0] = iD[0]
        D[1] = iD[1]
        D[2] = iD[2]
        D[3] = jD[0]
        D[4] = jD[1]
        D[5] = jD[2]

        return D

    def Dlocal(self, combo):
        """

        Parameters
        ----------
        combo : TYPE
            DESCRIPTION.

        Returns
        -------
        Dlocal : TYPE
            DESCRIPTION.

        """

        Dglobal = self.Dglobal(combo)

        Dlocal = np.matmul(self.T(), Dglobal)

        return Dlocal

    def Flocal(self, combo):
        """

        Parameters
        ----------
        combo : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        Dlocal = self.Dlocal(combo)
        Qf = np.reshape(self.FEF(combo), (-1, 1))

        FL = np.matmul(self.k(), Dlocal.T)

        self.end_forces_local[combo] = FL + Qf

    def Fglobal(self, combo):
        """

        Parameters
        ----------
        combo : TYPE
            DESCRIPTION.

        Returns
        -------
        FG : TYPE
            DESCRIPTION.

        """

        Dglobal = self.Dglobal(combo)
        Qfg = self.FEFglobal(combo)

        # global stiffness matrix
        KG = self.kglobal()

        FG = np.matmul(KG, Dglobal)

        self.end_forces_global[combo] = FG + Qfg

        self.Flocal(combo)

        return FG + Qfg

    def stations(self, num_stations=10):

        """
        define general computation points along the beam length for shear,
        moment, slope, and deflection plots
        """

        # parametric list of stations between 0 and 1'
        eta = [0 + i * (1 / num_stations) for i in range(num_stations + 1)]

        stations = [self.length * i for i in eta]

        if self._loaded:
            extra_stations = []

            for load in self.loads:
                if (
                    load.kind == "POINT"
                    or load.kind == "MOMENT"
                    or load.kind == "AXIAL_POINT"
                ):
                    b = min(self.length, load.a + 0.001)
                    c = max(0, load.a - 0.001)
                    extra_stations.extend([c, load.a, b])

                elif load.kind == "LINE" or load.kind == "AXIAL_LINE":
                    c = min(self.length, load.b + 0.001)
                    d = max(0, load.a - 0.001)
                    extra_stations.extend([d, load.a, load.b, c])
                else:
                    pass

            stations.extend(extra_stations)

        stations.sort()

        # Make sure the first and last stations do not exceed the beam

        if stations[0] < 0:
            stations[0] = 0

        if stations[-1] > self.length:
            stations[-1] = self.length

        # Remove duplicate locations
        self.calcstations = sorted(set(stations))

        self._stations = True

    def dlocal_plot(self, combo, scale=1):

        if not self._stations:
            self.stations()

        Dlocal = self.Dlocal(combo)

        # Parametric Functions defining a linear relationship for deflection
        # in each axis based on the Ux and Uy nodal displacements
        Dx = lambda x: Dlocal[0, 0] + (x / self.length) * (
            Dlocal[0, 3] - Dlocal[0, 0]
        )
        Dy = lambda x: Dlocal[0, 1] + (x / self.length) * (
            Dlocal[0, 4] - Dlocal[0, 1]
        )

        empty_f = np.zeros((6, 1))

        Fendlocal = self.end_forces_local.get(combo, empty_f)

        # Empty Piecwise functions to build the total function from the loading
        dx = loadtypes.Piecewise_Polynomial()
        dy = loadtypes.Piecewise_Polynomial()

        # Create "loads" from the end forces and combine with dx and dy
        fxi = loadtypes.R2_Axial_Load(Fendlocal[0, 0], 0, self)
        fyi = loadtypes.R2_Point_Load(Fendlocal[1, 0], 0, self)
        mzi = loadtypes.R2_Point_Moment(Fendlocal[2, 0], 0, self)
        fxj = loadtypes.R2_Axial_Load(Fendlocal[3, 0], self.length, self)
        fyj = loadtypes.R2_Point_Load(Fendlocal[4, 0], self.length, self)
        mzj = loadtypes.R2_Point_Moment(Fendlocal[5, 0], self.length, self)

        dx = dx.combine(fxi.Dx, 1, 1)
        dy = dy.combine(fyi.Dy, 1, 1)
        dy = dy.combine(mzi.Dy, 1, 1)
        dx = dx.combine(fxj.Dx, 1, 1)
        dy = dy.combine(fyj.Dy, 1, 1)
        dy = dy.combine(mzj.Dy, 1, 1)

        # Combine Piecewise Deflection Functions of all of the loads
        if self._loaded:

            for load in self.loads:

                if load.loadcase == combo:

                    dx = dx.combine(load.Dx, 1, 1)
                    dy = dy.combine(load.Dy, 1, 1)

        dlocal_span = np.zeros((len(self.calcstations), 2))

        for i, x in enumerate(self.calcstations):

            dxl = dx.evaluate(x) + Dx(x)
            dyl = dy.evaluate(x) + Dy(x)

            dlocal_span[i, 0] = x + (dxl * scale)
            dlocal_span[i, 1] = dyl * scale

        return dlocal_span

    def dglobal_plot(self, combo, scale=1):

        dlocal_plot = self.dlocal_plot(combo, scale=scale)

        c = (self.jnode.x - self.inode.x) / self.length
        s = (self.jnode.y - self.inode.y) / self.length

        R = np.matrix([[c, s], [-s, c]])

        dglobal_plot = np.matmul(dlocal_plot, R)

        return dglobal_plot

    def Mlocal_plot(self, combo, scale=1):

        if not self._stations:
            self.stations()

        empty_f = np.zeros((6, 1))

        Fendlocal = self.end_forces_local.get(combo, empty_f)

        # Empty Piecwise functions to build the total function from the loading
        Mzx = loadtypes.Piecewise_Polynomial()

        # Create "loads" from the end forces and combine with dx and dy
        fxi = loadtypes.R2_Axial_Load(Fendlocal[0, 0], 0, self)
        fyi = loadtypes.R2_Point_Load(Fendlocal[1, 0], 0, self)
        mzi = loadtypes.R2_Point_Moment(Fendlocal[2, 0], 0, self)
        fxj = loadtypes.R2_Axial_Load(Fendlocal[3, 0], self.length, self)
        fyj = loadtypes.R2_Point_Load(Fendlocal[4, 0], self.length, self)
        mzj = loadtypes.R2_Point_Moment(Fendlocal[5, 0], self.length, self)

        Mzx = Mzx.combine(fxi.Mz, 1, 1)
        Mzx = Mzx.combine(fyi.Mz, 1, 1)
        Mzx = Mzx.combine(mzi.Mz, 1, 1)
        Mzx = Mzx.combine(fxj.Mz, 1, 1)
        Mzx = Mzx.combine(fyj.Mz, 1, 1)
        Mzx = Mzx.combine(mzj.Mz, 1, 1)

        # Combine Piecewise Deflection Functions of all of the loads
        if self._loaded:

            for load in self.loads:

                if load.loadcase == combo:

                    Mzx = Mzx.combine(load.Mz, 1, 1)

        mlocal_span = np.zeros((len(self.calcstations), 2))

        for i, x in enumerate(self.calcstations):

            m = Mzx.evaluate(x)

            mlocal_span[i, 0] = x
            mlocal_span[i, 1] = m * scale

        return mlocal_span

    def Mglobal_plot(self, combo, scale):

        mlocal_plot = self.Mlocal_plot(combo, scale=scale)

        c = (self.jnode.x - self.inode.x) / self.length
        s = (self.jnode.y - self.inode.y) / self.length

        R = np.matrix([[c, s], [-s, c]])

        mglobal_plot = np.matmul(mlocal_plot, R)

        return mglobal_plot
