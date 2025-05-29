# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 9), (2, 3)]
colors = [1, 2, 3]

# Initialize the coloring map
coloring = {}

# Assign colors to the connected vertices
coloring[1] = 1
coloring[9] = 2
coloring[2] = 1
coloring[3] = 2

# Assign colors to the isolated vertices
for vertex in vertices:
    if vertex not in coloring:
        coloring[vertex] = 1

# Print the coloring map
print(coloring)