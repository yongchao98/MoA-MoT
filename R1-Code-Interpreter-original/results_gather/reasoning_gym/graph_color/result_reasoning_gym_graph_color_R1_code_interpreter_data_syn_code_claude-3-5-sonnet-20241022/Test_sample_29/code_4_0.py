import json

# Define the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (0, 8), (1, 4), (1, 6), (1, 8), (5, 6), (5, 8)]

# Create a solution by analyzing the graph structure
# Looking at the edges, we can see that:
# - Vertices 6 and 8 are the most connected (hub vertices)
# - Vertices 0, 1, and 5 are connected to these hubs
# - The rest of the vertices have no connections

def create_valid_coloring():
    coloring = {}
    
    # First, assign color 3 to one of the hub vertices (6)
    coloring["6"] = 3
    
    # Assign color 2 to the other hub vertex (8)
    coloring["8"] = 2
    
    # Vertices connected to hubs (0, 1, 5) must have color 1
    coloring["0"] = 1
    coloring["1"] = 1
    coloring["5"] = 1
    
    # Vertex 4 is connected to vertex 1, so it needs a different color
    coloring["4"] = 2
    
    # Remaining vertices (2, 3, 7, 9) have no connections, can be color 1
    coloring["2"] = 1
    coloring["3"] = 1
    coloring["7"] = 1
    coloring["9"] = 1
    
    return coloring

# Create and verify the solution
coloring = create_valid_coloring()

# Verify that no adjacent vertices have the same color
for edge in edges:
    v1, v2 = str(edge[0]), str(edge[1])
    if coloring[v1] == coloring[v2]:
        raise Exception(f"Invalid coloring: vertices {v1} and {v2} have same color {coloring[v1]}")

print(f"<<<{json.dumps(coloring)}>>>")