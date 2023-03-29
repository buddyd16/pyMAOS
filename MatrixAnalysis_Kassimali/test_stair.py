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

# Units
# Length - inches
# Force - kips

# Sloped beam to test braced frame
loadcase = "D"
loadcase2 = "L"
loadcombo2 = LoadCombo("S1", {"D": 1}, ["D"], False, "SLS")
loadcombo = LoadCombo("S2", {"L": 1}, ["L"], False, "SLS")

# Nodes
N1 = R2Node(0, 0)
N2 = R2Node(6.16 * 12, 5.76 * 12)
N3 = R2Node((6.16 + 3.67) * 12, 5.76 * 12)
N4 = R2Node((6.16 + 3.67 + 6.33) * 12, (5.76 + 5.1) * 12)


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
CIPMaterial = Material((110 / (12 * 12 * 12.0)), 2408)


# Sections
# 9x44
StairSection = Section(396, 2673, 63888)


# Members
RF1 = R2Frame(N1, N2, CIPMaterial, StairSection, 1)
RF2 = R2Frame(N2, N3, CIPMaterial, StairSection, 2)
RF3 = R2Frame(N3, N4, CIPMaterial, StairSection, 3)

# Member List
members = [RF1, RF2, RF3]

# Member Release
# RF2.hinge_i()
# RF2.hinge_j()

# Member Loads
sw = -0.3025 / 12
RF1.add_distributed_load(sw, sw, 0, 100, loadcase, "Y", location_percent=True)
RF2.add_distributed_load(sw, sw, 0, 100, loadcase, "Y", location_percent=True)
RF3.add_distributed_load(sw, sw, 0, 100, loadcase, "Y", location_percent=True)

sdl = -0.10093 / 12
RF1.add_distributed_load(
    sdl, sdl, 0, 100, loadcase, "Y", location_percent=True
)
RF2.add_distributed_load(
    sdl, sdl, 0, 100, loadcase, "Y", location_percent=True
)
RF3.add_distributed_load(
    sdl, sdl, 0, 100, loadcase, "Y", location_percent=True
)

ll = (-100 * 3.67) / (12 * 1000)
RF1.add_distributed_load(
    ll, ll, 0, 100, loadcase2, "Y", location_percent=True, projected=True
)
RF2.add_distributed_load(
    ll, ll, 0, 100, loadcase2, "Y", location_percent=True, projected=True
)
RF3.add_distributed_load(
    ll, ll, 0, 100, loadcase2, "Y", location_percent=True, projected=True
)

# Create the 2D Structure
Structure = R2Struct.R2Structure(nodes, members)


Structure.set_node_uids()
Structure.set_member_uids()
FM = Structure.freedom_map()
K = Structure.Kstructure(FM)
U = Structure.solve_linear_static(loadcombo)
U2 = Structure.solve_linear_static(loadcombo2)
Errors = Structure._ERRORS

# Print Output
print("Errors:")
print(Errors)
print(f"Loadcase: {loadcase}")
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
print("-" * 100)
print(f"Loadcase: {loadcase2}")
print("Displacements:")
for i, node in enumerate(nodes):
    tx = node.displacements[loadcombo2.name]
    print(
        f"N{node.uid} -- Ux: {tx[0]:.4E} -- Uy:{tx[1]:.4E} -- Rz:{tx[2]:.4E}"
    )
print("-" * 100)
print("Reactions:")
for i, node in enumerate(nodes):
    rx = node.reactions[loadcombo2.name]
    print(
        f"N{node.uid} -- Rx: {rx[0]:.4E} -- Ry:{rx[1]:.4E} -- Mz:{rx[2]:.4E}"
    )
print("-" * 100)
print("Member Forces:")
for i, member in enumerate(members):
    fx = member.end_forces_local[loadcombo2.name]

    print(f"M{member.uid}")
    print(
        f"    i -- Axial: {fx[0,0]:.4E} -- Shear: {fx[1,0]:.4E} -- Moment: {fx[2,0]:.4E}"
    )
    print(
        f"    j -- Axial: {fx[3,0]:.4E} -- Shear: {fx[4,0]:.4E} -- Moment: {fx[5,0]:.4E}"
    )

# Plot the structure
fig, axs = plt.subplots(2, 3, figsize=(8, 8))

axial_scale = 1
shear_scale = 10
moment_scale = 0.2
rotation_scale = 10000
displace_scale = 100

axs[0, 0].set_title(
    f"Geometry and Deformed Shape\n scale:{displace_scale}", fontsize=12
)

axs[0, 1].set_title(f"Axial Force\n scale:{axial_scale}", fontsize=12)

axs[0, 2].set_title(f"Shear Force\n scale:{shear_scale}", fontsize=12)

axs[1, 0].set_title(f"Moment\n scale:{moment_scale}", fontsize=12)

axs[1, 1].set_title(
    f"Cross-Section Rotation\n scale:{rotation_scale}", fontsize=12
)

for node in nodes:
    axs[0, 0].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[0, 1].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[0, 2].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[1, 0].plot(node.x, node.y, marker=".", markersize=8, color="red")
    axs[1, 1].plot(node.x, node.y, marker=".", markersize=8, color="red")
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
    axs[0, 2].plot(
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

    axs[0, 2].plot(
        (vglobal[:, 0] + member.inode.x),
        (vglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="green",
    )

    axs[1, 0].plot(
        (mglobal[:, 0] + member.inode.x),
        (mglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="red",
    )
    axs[1, 0].plot(
        (member.inode.x, mglobal[0, 0] + member.inode.x),
        (member.inode.y, mglobal[0, 1] + member.inode.y),
        linewidth=1,
        color="red",
    )
    axs[1, 0].plot(
        (member.jnode.x, mglobal[-1, 0] + member.inode.x),
        (member.jnode.y, mglobal[-1, 1] + member.inode.y),
        linewidth=1,
        color="red",
    )

    axs[0, 0].plot(
        (dglobal[:, 0] + member.inode.x),
        (dglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="gray",
    )

    axs[1, 1].plot(
        (sglobal[:, 0] + member.inode.x),
        (sglobal[:, 1] + member.inode.y),
        linewidth=1,
        color="purple",
    )

axs[0, 0].grid(True)
axs[0, 1].grid(True)
axs[0, 2].grid(True)
axs[1, 0].grid(True)
axs[1, 1].grid(True)
axs[1, 2].grid(True)

axs[0, 0].set_aspect("equal", "box")
axs[0, 1].set_aspect("equal", "box")
axs[0, 2].set_aspect("equal", "box")
axs[1, 0].set_aspect("equal", "box")
axs[1, 1].set_aspect("equal", "box")
axs[1, 2].set_aspect("equal", "box")

fig.tight_layout()

plt.show()
