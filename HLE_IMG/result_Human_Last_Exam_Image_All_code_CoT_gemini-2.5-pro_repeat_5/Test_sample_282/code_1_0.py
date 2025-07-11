# Define the directed graph as an adjacency list
# The keys are the vertices, and the values are lists of vertices they point to.
graph = {
    0: [3],
    1: [8],
    2: [4],
    3: [2, 6],
    4: [1, 5, 6],
    5: [2, 8],
    6: [1, 5],
    7: [],
    8: [7]
}

# The sum of out-degrees is the total number of edges in the graph.
# We can calculate it by summing the lengths of the adjacency lists for each vertex.
total_edges = 0
calculation_steps = []

for vertex, neighbors in sorted(graph.items()):
    out_degree = len(neighbors)
    total_edges += out_degree
    calculation_steps.append(str(out_degree))

# The term "degree sum" in statement E likely refers to the number of edges,
# which is the sum of all out-degrees.
print("Calculating the total number of edges (sum of out-degrees):")
print(" + ".join(calculation_steps), f"= {total_edges}")
print(f"The total number of edges is {total_edges}.")
print("Assuming 'degree sum' refers to the number of edges, this matches statement E.")
