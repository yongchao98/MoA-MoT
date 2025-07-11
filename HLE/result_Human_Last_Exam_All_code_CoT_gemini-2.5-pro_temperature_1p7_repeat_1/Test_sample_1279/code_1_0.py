import math

# --- Introduction ---
print("Calculating the symmetry factors for second-order vacuum bubble diagrams in φ⁴ theory.")
print("There are two topologically distinct connected diagrams at this order.\n")


# --- Diagram 1: The "Sunset" Diagram ---
print("="*40)
print("Diagram 1: The 'Sunset' Diagram")
print("="*40)
print("Description: Two vertices connected to each other by four propagators.")
print("\nSources of symmetry:")
print("1. Swapping the two identical vertices. (Factor: 2)")
print("2. Permuting the four identical propagators between the vertices. (Factor: 4!)")

# Calculation for Diagram 1
vert_sym_1 = 2
prop_sym_1 = math.factorial(4)
s1 = vert_sym_1 * prop_sym_1

print("\nSymmetry factor calculation (S1):")
print(f"S1 = (vertex symmetry) * (propagator symmetry)")
print(f"S1 = {vert_sym_1} * {prop_sym_1} = {s1}")
print("\n")


# --- Diagram 2: The "Double-Scoop" Diagram ---
print("="*40)
print("Diagram 2: The 'Double-Scoop' Diagram")
print("="*40)
print("Description: Each of the two vertices has one self-loop. The vertices are")
print("             connected to each other by two propagators.")
print("\nSources of symmetry:")
print("1. Swapping the two identical vertices (with their loops). (Factor: 2)")
print("2. Swapping the two identical propagators that connect the vertices. (Factor: 2)")
print("3. Swapping the ends of the propagator forming the first self-loop. (Factor: 2)")
print("4. Swapping the ends of the propagator forming the second self-loop. (Factor: 2)")

# Calculation for Diagram 2
vert_sym_2 = 2
conn_prop_sym_2 = 2
loop1_sym_2 = 2
loop2_sym_2 = 2
s2 = vert_sym_2 * conn_prop_sym_2 * loop1_sym_2 * loop2_sym_2

print("\nSymmetry factor calculation (S2):")
print(f"S2 = (vertex sym) * (connecting propagator sym) * (loop 1 sym) * (loop 2 sym)")
print(f"S2 = {vert_sym_2} * {conn_prop_sym_2} * {loop1_sym_2} * {loop2_sym_2} = {s2}")
print("\n")


# --- Total Sum ---
print("="*40)
print("Total Sum of Symmetry Factors")
print("="*40)
total_sum = s1 + s2
print("The total is the sum of the symmetry factors for each diagram.")
print(f"Total Sum = S1 + S2")
print(f"Total Sum = {s1} + {s2} = {total_sum}")

<<<64>>>