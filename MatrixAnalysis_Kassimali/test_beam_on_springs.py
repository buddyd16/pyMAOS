from nodes import R2Node
from elements import R2Truss, R2Frame
from material import LinearElasticMaterial as Material
from section import Section
import R2Structure as R2Struct
from loadcombos import LoadCombo
import plot_structure


########################################
#                                      #
# CONSISTENT UNIT SYSTEM REQUIRED !!!! #
#                                      #
########################################

# Units
# Length - inches
# Force - kips

loadcase = "D"
loadcombo = LoadCombo("S1", {"D": 1}, ["D"], False, "SLS")

# Nodes
num_nodes = 41
# parametric list of between 0 and 1'
eta = [0 + i * (1 / num_nodes) for i in range(num_nodes)]
length_ft = 40
length_in = length_ft * 12
# Generate Nodes
nodes = []
for pt in eta:
    x = pt*length_in
    nodes.append(R2Node(x, 0))

# Apply the springs
k = 10 #kpi

for node in nodes:
    node.releaseAll()
    node.applySpringUy(k,0)

# Restrain last node for Ux to prevent rigid body motion in sliding
nodes[-1].restrainUx()

# Nodal Loads
position = 19
nodes[position].loads[loadcase] = [0, -1000, 0]

# Materials
BeamMaterial = Material(0.00028, 29000)

# Sections
# W12x14
BeamSection = Section(4.16, 88.6,2.36)

members = []
# Generate Members
for i, node in enumerate(nodes):
    if i == 0:
        pass
    else:
        members.append(R2Frame(nodes[i-1],node,BeamMaterial,BeamSection))

# Create the 2D Structure
Structure = R2Struct.R2Structure(nodes, members)
Structure.set_node_uids()
Structure.set_member_uids()

Structure.spring_nodes()

FM = Structure.freedom_map()
K = Structure.Kstructure(FM)
U = Structure.solve_linear_static(loadcombo)
Errors = Structure._ERRORS
print(Errors)

# Plot the structure
scaling = {
        "axial_load": 1,
        "normal_load": 1,
        "point_load": 1,
        "axial": 1,
        "shear": 1,
        "moment": 0.001,
        "rotation": 100,
        "displacement": 10,
    }

plot_structure.plot_structure(nodes, members, loadcombo, scaling)