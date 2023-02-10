from nodes import R2Node
from elements import R2Truss, R2Frame
from material import LinearElasticMaterial as Material
from section import Section
import loading as loads
import R2Structure as R2Struct

import matplotlib.pyplot as plt


##########################################
##                                      ##
## CONSISTENT UNIT SYTEM REQUIRED !!!!! ##
##                                      ##
##########################################


# Material and Section
M1 = Material(29000)
M2 = Material(10000)
S1 = Section(8, 2370)
S2 = Section(12, 2370)
S3 = Section(16, 2370)

# Nodes
N1 = R2Node(0, 0, 1)
N2 = R2Node(288, 0, 2)
N3 = R2Node(576, 0, 3)
N4 = R2Node(864, 0, 4)
N5 = R2Node(288, 216, 5)
N6 = R2Node(576, 216, 6)

# Restraints
# N1.restraints = [1,1,1]
# N2.restraints = [0,0,1]
# N3.restraints = [0,1,1]
# N4.restraints = [0,1,1]
# N5.restraints = [0,0,1]
# N6.restraints = [0,0,1]

N1.restraints = [1, 1, 0]
N2.restraints = [0, 0, 0]
N3.restraints = [0, 1, 0]
N4.restraints = [0, 1, 0]
N5.restraints = [0, 0, 0]
N6.restraints = [0, 0, 0]

# Loading
loadcase = "D"

N2.loads[loadcase] = [0, -75, 0]
N5.loads[loadcase] = [25, 0, 0]
N6.loads[loadcase] = [0, -60, 0]

# Node List
nodes = [N1, N2, N3, N4, N5, N6]


# Elements
# T1 = R2Truss(N1, N2, M1, S1, 1)
# T2 = R2Truss(N2, N3, M1, S1, 2)
# T3 = R2Truss(N3, N4, M2, S3, 3)
# T4 = R2Truss(N5, N6, M1, S1, 4)
# T5 = R2Truss(N2, N5, M1, S1, 5)
# T6 = R2Truss(N3, N6, M1, S1, 6)
# T7 = R2Truss(N1, N5, M1, S2, 7)
# T8 = R2Truss(N2, N6, M1, S2, 8)
# T9 = R2Truss(N3, N5, M1, S2, 9)
# T10 = R2Truss(N4, N6, M2, S3, 10)

T1 = R2Frame(N1, N2, M1, S1, 1)
T2 = R2Frame(N2, N3, M1, S1, 2)
T3 = R2Frame(N3, N4, M2, S3, 3)
T4 = R2Frame(N5, N6, M1, S1, 4)
T5 = R2Frame(N2, N5, M1, S1, 5)
T6 = R2Frame(N3, N6, M1, S1, 6)
T7 = R2Frame(N1, N5, M1, S2, 7)
T8 = R2Frame(N2, N6, M1, S2, 8)
T9 = R2Frame(N3, N5, M1, S2, 9)
T10 = R2Frame(N4, N6, M2, S3, 10)

# Member List
members = [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10]

# Create the 2D Structure
Structure = R2Struct.R2Structure(nodes, members)

Errors = Structure._ERRORS
K = Structure.Kstructure()
U = Structure.solve_linear_static("D")

print(Errors)
print("Displacements:")
for node in nodes:
    tx = node.displacements["D"]
    print(f"N{node.uid} -- Ux: {tx[0]:.4E}  Uy:{tx[1]:.4E}  Rz:{tx[2]:.4E}")
print("-" * 100)
print("Reactions:")
for node in nodes:
    rx = node.reactions["D"]
    print(f"N{node.uid} -- Rx: {rx[0]:.4E}  Ry:{rx[1]:.4E}  Mz:{rx[2]:.4E}")
print("-" * 100)
print("Member Forces:")
for member in members:
    fx = member.end_forces_local["D"]
    print(f"M{member.uid} -- Axial: {fx[0,0]:.4E}")
# plot the structure
fig, ax = plt.subplots()

for node in nodes:
    ax.plot(node.x, node.y, marker=".", markersize=20, color="red")
    ax.plot(
        node.x_displaced("D", 100),
        node.y_displaced("D", 100),
        marker=".",
        markersize=10,
        color="gray",
    )
for member in members:
    ax.plot(
        [member.inode.x, member.jnode.x],
        [member.inode.y, member.jnode.y],
        linewidth=1,
        color="blue",
    )
    ax.plot(
        [member.inode.x_displaced("D", 100), member.jnode.x_displaced("D", 100)],
        [member.inode.y_displaced("D", 100), member.jnode.y_displaced("D", 100)],
        linewidth=1,
        color="gray",
    )
ax.grid(True)
fig.tight_layout()

plt.show()

# Load test bed
# N1 = R2Node(0, 0, 1)
# N2 = R2Node(10, 0, 2)

# M5 = Material(4176000.0)
# S5 = Section(0.02056,0.0014853395061728396)

# B1 = R2Frame(N1, N2, M5, S5, 1)

# lineload = loads.R2_Axial_Linear_Load(-10, -5, 2, 8, B1)
# ptload = loads.R2_Axial_Load(24, 0, B1)
# ptload2 = loads.R2_Axial_Load(21, B1.length, B1)

# print(lineload.Rix)
# print(lineload.Rjx)
# print(ptload.Rix)
# print(ptload.Rjx)
# print(ptload2.Rix)
# print(ptload2.Rjx)
# print(lineload.Ax)
# print(lineload.Dx)
# print(ptload.Ax)
# print(ptload.Dx)
# print(ptload2.Ax)
# print(ptload2.Dx)
# print(lineload.FEF())
# print(ptload.FEF())
# print(ptload2.FEF())


# A = ptload2.Ax.combine(ptload.Ax.combine(lineload.Ax, 1, 1),1,1)
# D = ptload2.Dx.combine(ptload.Dx.combine(lineload.Dx, 12, 12),12,1)

# num_stations = 100
# eta = [0+i*(1/num_stations)*B1.length for i in range(num_stations+1)]
# v = [A.evaluate(i) for i in eta]
# fig, ax = plt.subplots()
# ax.plot(eta,v, marker=".", markersize=2, color='red')
# ax.grid(True)
# fig.tight_layout()

# plt.show()


N7 = R2Node(0, 0, 7)
N8 = R2Node(120, 240, 8)
N9 = R2Node(360, 240, 9)

M4 = Material(29000)
S4 = Section(11.8, 310)

M11 = R2Frame(N7, N8, M4, S4, 11)
M12 = R2Frame(N8, N9, M4, S4, 12)

M11.add_point_load(-90, 50, "D", "Y", location_percent=True)
M12.add_distributed_load(-0.125, -0.125, 0, 100, "D", "Y", location_percent=True)

for load in M11.loads:
    print(load.p)
FEF = M11.FEF("D")
FEFg = M11.FEFglobal("D")

FEFm12 = M12.FEF("D")
FEFgm12 = M12.FEFglobal("D")

print(FEFg[0, 3:])
print(FEFgm12[0, 0:3])
Pf = FEFg[0, 3:] + FEFgm12[0, 0:3]
FEFg = M11.FEFglobal("D")

FEFm12 = M12.FEF("D")
FEFgm12 = M12.FEFglobal("D")

print(FEFg[0, 3:])
print(FEFgm12[0, 0:3])
Pf = FEFg[0, 3:] + FEFgm12[0, 0:3]
