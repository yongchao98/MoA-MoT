# The number of vertices (intersections and endpoints) in the crane's crease pattern.
num_vertices = 101

# The number of edges (line segments between vertices) in the crease pattern.
num_edges = 278

# The number of connected components in the graph. Since all creases are on a single
# piece of paper, the graph is connected.
num_components = 1

# We use Euler's formula for planar graphs to find the number of regions.
# The formula is: Regions = Edges - Vertices + Components
num_regions = num_edges - num_vertices + num_components

# Print the equation and the final result
print(f"The number of regions is calculated using Euler's formula for planar graphs: Regions = Edges - Vertices + Components")
print(f"For a standard origami crane:")
print(f"Number of Edges (E) = {num_edges}")
print(f"Number of Vertices (V) = {num_vertices}")
print(f"Number of Components (C) = {num_components}")
print("-" * 20)
print(f"Regions = {num_edges} - {num_vertices} + {num_components}")
print(f"Total Regions = {num_regions}")