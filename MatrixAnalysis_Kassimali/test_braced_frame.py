from nodes import R2Node
from elements import R2Truss, R2Frame
from material import LinearElasticMaterial as Material
from section import Section
import R2Structure as R2Struct
from loadcombos import LoadCombo

import matplotlib.pyplot as plt


##########################################
##                                      ##
## CONSISTENT UNIT SYTEM REQUIRED !!!!! ##
##                                      ##
##########################################

# Sloped beam to test braced frame
loadcase = "D"
loadcombo = LoadCombo("S1", {"D": 1}, ["D"], False, "SLS")

# Nodes
N1 = R2Node(0, 0)
N2 = R2Node(0, 120)
N3 = R2Node(120, 120)
N4 = R2Node(120, 0)


# Node Restraints
N1.restraints = [1, 1, 0]
N2.restraints = [0, 0, 0]
N3.restraints = [0, 0, 0]
N4.restraints = [1, 1, 0]

# Node List
nodes = [N1, N2, N3, N4]

# Nodal Loads
# N2.loads[loadcase] = [50, 0, 0]


# Materials
SteelMaterial = Material(0.00028, 29000)


# Sections
# W16x40
BeamSection = Section(11.8, 518, 28.9)
# W8x24, weak
ColumnSection = Section(7.08, 18.3, 82.7)
# L2x2x3/8
BraceSection = Section(1.37, 0.476, 0.476)

# Members
RF1 = R2Frame(N1, N2, SteelMaterial, ColumnSection, 1)
RF2 = R2Frame(N2, N3, SteelMaterial, BeamSection, 2)
RF3 = R2Frame(N3, N4, SteelMaterial, ColumnSection, 3)
RF4 = R2Truss(N2, N4, SteelMaterial, BraceSection, 4)
RF5 = R2Truss(N1, N3, SteelMaterial, BraceSection, 5)

# Member List
members = [RF1, RF2, RF3, RF4, RF5]
# members = [RF1, RF2, RF3, RF5]

# Member Release
RF2.hinge_i()
RF2.hinge_j()

# Member Loads
RF1.add_distributed_load(
    1 / 12, 1 / 12, 0, 100, loadcase, "X", location_percent=True
)
RF2.add_distributed_load(
    -1 / 12, -1 / 12, 0, 100, loadcase, "yy", location_percent=True
)
RF3.add_distributed_load(
    0.5 / 12, 0.5 / 12, 0, 100, loadcase, "X", location_percent=True
)

# Create the 2D Structure
Structure = R2Struct.R2Structure(nodes, members)


Structure.set_node_uids()
Structure.set_member_uids()
FM = Structure.freedom_map()
K = Structure.Kstructure(FM)
U = Structure.solve_linear_static(loadcombo)
Errors = Structure._ERRORS

# Print Output
print("Errors:")
print(Errors)
print("Displacements:")
for i, node in enumerate(nodes):
    tx = node.displacements[loadcombo.name]
    print(
        f"N{node.uid} -- Ux: {tx[0]:.4E} -- Uy:{tx[1]:.4E} -- Rz:{tx[2]:.4E}"
    )
print("-" * 100)
print("Reactions:")
for i, node in enumerate(nodes):
    rx = node.reactions[loadcombo.name]
    print(
        f"N{node.uid} -- Rx: {rx[0]:.4E} -- Ry:{rx[1]:.4E} -- Mz:{rx[2]:.4E}"
    )
print("-" * 100)
print("Member Forces:")
for i, member in enumerate(members):
    fx = member.end_forces_local[loadcombo.name]

    print(f"M{member.uid}")
    print(
        f"    i -- Axial: {fx[0,0]:.4E} -- Shear: {fx[1,0]:.4E} -- Moment: {fx[2,0]:.4E}"
    )
    print(
        f"    j -- Axial: {fx[3,0]:.4E} -- Shear: {fx[4,0]:.4E} -- Moment: {fx[5,0]:.4E}"
    )


# Plot the structure
fig, axs = plt.subplots(3, 2)

axial_scale = 0.5
shear_scale = 1
moment_scale = 0.25
rotation_scale = 500
displace_scale = 10

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
        node.x_displaced(loadcombo, displace_scale),
        node.y_displaced(loadcombo, displace_scale),
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
    aglobal = member.Aglobal_plot(loadcombo, axial_scale)
    vglobal = member.Vglobal_plot(loadcombo, shear_scale)
    mglobal = member.Mglobal_plot(loadcombo, moment_scale)
    sglobal = member.Sglobal_plot(loadcombo, rotation_scale)
    dglobal = member.Dglobal_plot(loadcombo, displace_scale)

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

plt.axis("square")

fig.tight_layout()

plt.show()
