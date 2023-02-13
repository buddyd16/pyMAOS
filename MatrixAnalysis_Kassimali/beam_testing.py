from nodes import R2Node
from elements import R2Truss, R2Frame
from material import LinearElasticMaterial as Material
from section import Section
import R2Structure as R2Struct

import matplotlib.pyplot as plt


##########################################
##                                      ##
## CONSISTENT UNIT SYTEM REQUIRED !!!!! ##
##                                      ##
##########################################

# 5 Element Cantilever beam w/ Point Load
loadcase = "D"

# Nodes
N1 = R2Node(0, 0, 1)
N2 = R2Node(60, 0, 2)
N3 = R2Node(120, 0, 3)
N4 = R2Node(180, 0, 4)
N5 = R2Node(240, 0, 5)
N6 = R2Node(300, 0, 6)

N7 = R2Node(0, 30, 7)
N8 = R2Node(300, 30, 8)

# Node Restraints
N1.restraints = [1, 1, 1]
N7.restraints = [1, 1, 1]

# Node List
nodes = [N1, N2, N3, N4, N5, N6, N7, N8]

# Nodal Loads
N6.loads[loadcase] = [0, -10, 0]
N8.loads[loadcase] = [0, -10, 0]

# Materials
BeamMaterial = Material(29000)


# Sections
# W24x55
BeamSection = Section(16.2, 1350)

# Members
RF1 = R2Frame(N1, N2, BeamMaterial, BeamSection, 1)
RF2 = R2Frame(N2, N3, BeamMaterial, BeamSection, 2)
RF3 = R2Frame(N3, N4, BeamMaterial, BeamSection, 3)
RF4 = R2Frame(N4, N5, BeamMaterial, BeamSection, 4)
RF5 = R2Frame(N5, N6, BeamMaterial, BeamSection, 5)
RF6 = R2Frame(N7, N8, BeamMaterial, BeamSection, 6)

# Member List
members = [RF1, RF2, RF3, RF4, RF5, RF6]

# Member Loads

# Create the 2D Structure
Structure = R2Struct.R2Structure(nodes, members)

Errors = Structure._ERRORS
K = Structure.Kstructure()
U = Structure.solve_linear_static(loadcase)

# Print Output
print("Errors:")
print(Errors)
print("Displacements:")
for i, node in enumerate(nodes):
    tx = node.displacements[loadcase]
    print(
        f"N{node.uid} -- Ux: {tx[0]:.4E} -- Uy:{tx[1]:.4E} -- Rz:{tx[2]:.4E}"
    )
print("-" * 100)
print("Reactions:")
for i, node in enumerate(nodes):
    rx = node.reactions[loadcase]
    print(
        f"N{node.uid} -- Rx: {rx[0]:.4E} -- Ry:{rx[1]:.4E} -- Mz:{rx[2]:.4E}"
    )
print("-" * 100)
print("Member Forces:")
for i, member in enumerate(members):
    fx = member.end_forces_local[loadcase]

    print(f"M{member.uid}")
    print(
        f"    i -- Axial: {fx[0,0]:.4E} -- Shear: {fx[1,0]:.4E} -- Moment: {fx[2,0]:.4E}"
    )
    print(
        f"    j -- Axial: {fx[3,0]:.4E} -- Shear: {fx[4,0]:.4E} -- Moment: {fx[5,0]:.4E}"
    )


# Plot the structure
fig, ax = plt.subplots()

moment_scale = 0.005
displace_scale = 10

for node in nodes:
    ax.plot(node.x, node.y, marker=".", markersize=20, color="red")
    ax.plot(
        node.x_displaced(loadcase, displace_scale),
        node.y_displaced(loadcase, displace_scale),
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
    mglobal = member.Mglobal_plot(loadcase, moment_scale)
    dglobal = member.dglobal_plot(loadcase, displace_scale)

    ax.plot(
        (mglobal[:, 0] + member.inode.x),
        (mglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="red",
    )

    ax.plot(
        (dglobal[:, 0] + member.inode.x),
        (dglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="gray",
    )
ax.grid(True)
fig.tight_layout()

plt.show()
