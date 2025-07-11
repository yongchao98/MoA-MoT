import math

# Plan:
# 1. Define the properties of the diagrams we are interested in.
#    - Theory: phi^4
#    - Order: 2nd order in coupling constant (2 vertices)
#    - Type: Vacuum bubbles (no external lines)
#    - We will sum the symmetry factors for all unique *connected* diagrams.

# 2. Calculate the symmetry factor for the first connected diagram (Figure-8).
#    This diagram has two vertices connected by four propagators.
#    The symmetry factor S is the product of:
#    - The number of ways to permute the two identical vertices: 2! = 2
#    - The number of ways to permute the four identical propagators: 4! = 24
s1_vertex_perms = math.factorial(2)
s1_line_perms = math.factorial(4)
s1 = s1_vertex_perms * s1_line_perms

print("Diagram 1: The 'Figure-8' Diagram")
print(f"Symmetry from permuting vertices = {s1_vertex_perms}")
print(f"Symmetry from permuting lines = {s1_line_perms}")
print(f"Symmetry factor S1 = {s1_vertex_perms} * {s1_line_perms} = {s1}\n")


# 3. Calculate the symmetry factor for the second connected diagram (Double-Scoop).
#    This diagram has two vertices, each with one loop, connected by two propagators.
#    The symmetries are:
#    - Swapping the two identical vertices (factor of 2).
#    - Swapping the two propagators connecting the vertices (factor of 2).
#    - Swapping the two lines that form the loop at the first vertex (factor of 2).
#    - Swapping the two lines that form the loop at the second vertex (factor of 2).
#    The total symmetry group size (the symmetry factor) for this diagram is 16.
s2 = 16
print("Diagram 2: The 'Double-Scoop' Diagram")
print("This diagram consists of two vertices, each with a self-loop,")
print("and two propagators connecting the vertices.")
print(f"The symmetry factor S2 for this diagram is {s2}\n")

# 4. Sum the symmetry factors of all connected second-order diagrams.
total_symmetry_factor = s1 + s2
print("The sum of the symmetry factors for all connected second-order vacuum bubbles is:")
print(f"Sum = S1 + S2 = {s1} + {s2} = {total_symmetry_factor}")
