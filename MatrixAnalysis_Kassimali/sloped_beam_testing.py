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
N2 = R2Node(8.66, 5, 2)
N3 = R2Node(20, 0, 3)
N4 = R2Node(28.66, 5, 4)
N5 = R2Node(0, 10, 5)
N6 = R2Node(8.66, 15, 6)
N7 = R2Node(20, 10, 7)
N8 = R2Node(28.66, 15, 8)

# Node Restraints
N1.restraints = [1, 1, 0]
N2.restraints = [1, 1, 0]
N3.restraints = [1, 1, 0]
N4.restraints = [1, 1, 0]
N5.restraints = [1, 1, 0]
N6.restraints = [1, 1, 0]
N7.restraints = [1, 1, 0]
N8.restraints = [1, 1, 0]

# Node List
nodes = [N1, N2, N3, N4, N5, N6, N7, N8]

# Nodal Loads

# Materials
BeamMaterial = Material(29000 * 6894.76)


# Sections
# W24x55
BeamSection = Section(0.0035, 0.000005770875286)

# Members
RF1 = R2Frame(N1, N2, BeamMaterial, BeamSection, 1)
RF2 = R2Frame(N3, N4, BeamMaterial, BeamSection, 2)
RF3 = R2Frame(N5, N6, BeamMaterial, BeamSection, 1)
RF4 = R2Frame(N7, N8, BeamMaterial, BeamSection, 2)

# Member List
members = [RF1, RF2, RF3, RF4]

# Member Release

# Member Loads
RF1.add_distributed_load(
    -1, -1, 0, 100, direction="Y", location_percent=True, projected=False
)
RF2.add_distributed_load(
    -1, -1, 0, 100, direction="Y", location_percent=True, projected=True
)
RF3.add_distributed_load(
    1, 1, 0, 100, direction="X", location_percent=True, projected=False
)
RF4.add_distributed_load(
    1, 1, 0, 100, direction="X", location_percent=True, projected=True
)

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
fig, axs = plt.subplots(3, 2)

axial_scale = 1.0
shear_scale = 1.0
moment_scale = 0.75
rotation_scale = 100
displace_scale = 50

axs[0, 0].set_title(
    f"Geometry and Deformed Shape\n scale:{displace_scale}", fontsize=12
)

axs[0, 1].set_title(f"Axial Force\n scale:{axial_scale}", fontsize=12)

axs[1, 0].set_title(f"Shear Force\n scale:{shear_scale}", fontsize=12)

axs[1, 1].set_title(f"Moment\n scale:{moment_scale}", fontsize=12)

axs[2, 0].set_title(
    f"Cross-Section Rotation\n scale:{rotation_scale}", fontsize=12
)

for node in nodes:
    axs[0, 0].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[0, 1].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[1, 0].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[1, 1].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[2, 0].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[0, 0].plot(
        node.x_displaced(loadcase, displace_scale),
        node.y_displaced(loadcase, displace_scale),
        marker=".",
        markersize=10,
        color="gray",
    )
for member in members:
    axs[0, 0].plot(
        [member.inode.x, member.jnode.x],
        [member.inode.y, member.jnode.y],
        linewidth=1,
        color="blue",
    )
    axs[0, 1].plot(
        [member.inode.x, member.jnode.x],
        [member.inode.y, member.jnode.y],
        linewidth=1,
        color="blue",
    )
    axs[1, 0].plot(
        [member.inode.x, member.jnode.x],
        [member.inode.y, member.jnode.y],
        linewidth=1,
        color="blue",
    )
    axs[1, 1].plot(
        [member.inode.x, member.jnode.x],
        [member.inode.y, member.jnode.y],
        linewidth=1,
        color="blue",
    )
    axs[2, 0].plot(
        [member.inode.x, member.jnode.x],
        [member.inode.y, member.jnode.y],
        linewidth=1,
        color="blue",
    )
    aglobal = member.Aglobal_plot(loadcase, axial_scale)
    vglobal = member.Vglobal_plot(loadcase, shear_scale)
    mglobal = member.Mglobal_plot(loadcase, moment_scale)
    sglobal = member.Sglobal_plot(loadcase, rotation_scale)
    dglobal = member.Dglobal_plot(loadcase, displace_scale)

    axs[0, 1].plot(
        (aglobal[:, 0] + member.inode.x),
        (aglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="blue",
    )

    axs[1, 0].plot(
        (vglobal[:, 0] + member.inode.x),
        (vglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="green",
    )

    axs[1, 1].plot(
        (mglobal[:, 0] + member.inode.x),
        (mglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="red",
    )

    axs[0, 0].plot(
        (dglobal[:, 0] + member.inode.x),
        (dglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="gray",
    )

    axs[2, 0].plot(
        (sglobal[:, 0] + member.inode.x),
        (sglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="purple",
    )

axs[0, 0].grid(True)
axs[0, 1].grid(True)
axs[1, 0].grid(True)
axs[1, 1].grid(True)
axs[2, 0].grid(True)

fig.tight_layout()

plt.show()