import math

# Description of the problem and method
print("This script calculates the symmetry factor (S) for a specific Feynman diagram.")
print("The diagram features two vertices and four field contractions (propagators) connecting them.")
print("\nThe symmetry factor is calculated based on the diagram's automorphisms (self-symmetries):")
print("1. Symmetry from permuting identical vertices.")
print("2. Symmetry from permuting identical propagators between the same two vertices.")

# Parameters of the diagram
num_vertices = 2
num_propagators = 4

print(f"\nFor this diagram:")
print(f"- There are {num_vertices} identical vertices that can be interchanged.")
print(f"- There are {num_propagators} identical propagators connecting these two vertices.")

# Calculate individual symmetry contributions
vertex_symmetry_factor = math.factorial(num_vertices)
propagator_symmetry_factor = math.factorial(num_propagators)

# Calculate the total symmetry factor
total_symmetry_factor = vertex_symmetry_factor * propagator_symmetry_factor

# Display the final calculation step-by-step
print("\nThe calculation for the total symmetry factor S is the product of these factors:")
print(f"S = (Symmetry of vertices) * (Symmetry of propagators)")
# The final response prints out each number in the equation.
print(f"S = {num_vertices}! * {num_propagators}!")
print(f"S = {vertex_symmetry_factor} * {propagator_symmetry_factor}")
print(f"S = {total_symmetry_factor}")