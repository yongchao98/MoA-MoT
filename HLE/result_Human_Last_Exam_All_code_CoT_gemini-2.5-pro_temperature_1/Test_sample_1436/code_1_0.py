import math

# Define the properties of the Feynman diagram
num_vertices = 2
num_propagators = 4

# Explain the problem and the method
print("This script calculates the symmetry factor for a Feynman diagram with two vertices connected by four propagators.")
print("-" * 80)
print("The symmetry factor (S) is determined by the number of ways the diagram can be rearranged to look the same.")
print("This is calculated from the permutations of its identical components.\n")

# Step 1: Calculate the symmetry factor from interchanging propagators
perms_propagators = math.factorial(num_propagators)
print(f"1. Symmetry from interchanging propagators:")
print(f"   There are {num_propagators} identical propagators connecting the same two vertices.")
print(f"   The number of permutations is {num_propagators}!, which is {perms_propagators}.\n")

# Step 2: Calculate the symmetry factor from interchanging vertices
perms_vertices = math.factorial(num_vertices)
print(f"2. Symmetry from interchanging vertices:")
print(f"   There are {num_vertices} identical vertices that can be swapped.")
print(f"   The number of permutations is {num_vertices}!, which is {perms_vertices}.\n")

# Step 3: Calculate the total symmetry factor
symmetry_factor = perms_propagators * perms_vertices
print("The total symmetry factor is the product of these two values.")
print("Final Equation:")
print(f"{num_propagators}! * {num_vertices}! = {symmetry_factor}")