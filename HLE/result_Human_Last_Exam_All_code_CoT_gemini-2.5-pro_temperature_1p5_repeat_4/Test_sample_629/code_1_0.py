import itertools

def perm_to_str(p):
    """Converts a 1-indexed permutation tuple to a string."""
    return f"({p[0]},{p[1]},{p[2]})"

def perm_mult(p1, p2):
    """
    Computes permutation multiplication p1 * p2 (p2 is applied first).
    Permutations are 1-indexed tuples. For example, in a 3x3 grid,
    p=(2,3,1) means 1->2, 2->3, 3->1.
    """
    n = len(p1)
    res = [0] * n
    # For each element i from 1 to n:
    # 1. Find where p2 sends it: p2[i-1]
    # 2. Find where p1 sends that result: p1[p2[i-1] - 1]
    for i in range(n):
        res[i] = p1[p2[i]-1]
    return tuple(res)

def perm_inv(p):
    """
    Computes the inverse of a permutation.
    If p sends i to p[i-1], the inverse sends p[i-1] back to i.
    """
    n = len(p)
    res = [0] * n
    for i in range(n):
        res[p[i]-1] = i + 1
    return tuple(res)

# Step 1: Set up the permutations for a 3x3 grid
S3 = sorted(list(itertools.permutations([1, 2, 3])))

# Step 2: Define the cycle for the left-hand trefoil
# The knot type is determined by O_inv * X. For the L-trefoil, this is a 3-cycle.
# By convention, we can choose the cycle c1 = (2, 3, 1), meaning 1->2, 2->3, 3->1.
c1 = (2, 3, 1)

print(f"There are {len(S3)}! = {len(S3)} permutations of 3 elements.")
print("A grid diagram (O, X) represents a left-hand trefoil if O_inv * X = c1.")
print(f"The chosen 3-cycle for the L-trefoil is c1 = {perm_to_str(c1)}.")
print("This implies that for any valid diagram (O, X), we have X = O * c1.")

# Step 3: Generate all 6 L-trefoil diagrams
l_trefoil_diagrams = []
for O in S3:
    X = perm_mult(O, c1)
    l_trefoil_diagrams.append((O, X))

print("\nThe 6 minimal grid diagrams for the L-trefoil are:")
for i, d in enumerate(l_trefoil_diagrams):
    print(f"  D{i+1}: (O={perm_to_str(d[0])}, X={perm_to_str(d[1])})")

# Step 4: Classify diagrams under translation (cyclic shifts)
cyclic_shifts = [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
orbits = []
remaining_diagrams = l_trefoil_diagrams.copy()

while remaining_diagrams:
    start_diagram = remaining_diagrams[0]
    current_orbit = {start_diagram}
    queue = [start_diagram]

    while queue:
        d = queue.pop(0)
        O, X = d
        # Action of translation group (p,q) on (O,X) is (p*O*q_inv, p*X*q_inv)
        for p in cyclic_shifts:
            for q in cyclic_shifts:
                q_inv = perm_inv(q)
                new_O = perm_mult(p, perm_mult(O, q_inv))
                new_X = perm_mult(p, perm_mult(X, q_inv))
                new_d = (new_O, new_X)
                if new_d in l_trefoil_diagrams and new_d not in current_orbit:
                    current_orbit.add(new_d)
                    queue.append(new_d)

    orbits.append(current_orbit)
    remaining_diagrams = [d for d in remaining_diagrams if d not in current_orbit]

print(f"\nUnder translations (cyclic shifts), the {len(l_trefoil_diagrams)} diagrams fall into {len(orbits)} orbits.")
for i, orbit in enumerate(orbits):
    orbit_size = len(orbit)
    print(f"Orbit {i+1} has {orbit_size} diagrams:")
    for d in sorted(list(orbit)):
        d_index = l_trefoil_diagrams.index(d) + 1
        print(f"  D{d_index}: (O={perm_to_str(d[0])}, X={perm_to_str(d[1])})")

print(f"\nEquation: {len(l_trefoil_diagrams)} total diagrams / {len(orbits[0])} diagrams per orbit = {len(orbits)} unique diagrams up to translation.")

# Step 5: Classify under rotation/reflection
print("\nNow we consider the full symmetry group of the square, which includes reflections.")
print("Let's test if a reflection maps one orbit to the other.")
# Reflection across the horizontal midline sends grid point (i,j) to (i, 4-j)
# This transforms a diagram (O,X) to (R*O, R*X) where R is the reversing permutation.
R = (3, 2, 1)
print(f"We apply reflection using the permutation R = {perm_to_str(R)}.")

d_from_orbit1 = list(orbits[0])[0]
d1_idx = l_trefoil_diagrams.index(d_from_orbit1) + 1
print(f"\nTake a diagram from Orbit 1, e.g., D{d1_idx}: (O={perm_to_str(d_from_orbit1[0])}, X={perm_to_str(d_from_orbit1[1])})")

reflected_O = perm_mult(R, d_from_orbit1[0])
reflected_X = perm_mult(R, d_from_orbit1[1])
reflected_d = (reflected_O, reflected_X)
rd_idx = l_trefoil_diagrams.index(reflected_d) + 1

print(f"Applying reflection to D{d1_idx} yields the new diagram (O',X'):")
print(f"  O' = R * O = {perm_to_str(R)} * {perm_to_str(d_from_orbit1[0])} = {perm_to_str(reflected_O)}")
print(f"  X' = R * X = {perm_to_str(R)} * {perm_to_str(d_from_orbit1[1])} = {perm_to_str(reflected_X)}")
print(f"This is diagram D{rd_idx}.")

is_in_orbit2 = reflected_d in orbits[1]
if is_in_orbit2:
    print(f"\nDiagram D{rd_idx} is in Orbit 2.")
    print("Since reflection connects Orbit 1 and Orbit 2, all 6 diagrams are equivalent under the full symmetry group.")
    print("\nFinal Conclusion: There is only 1 unique grid diagram up to translation and rotation.")
else:
    print("\nError in logic: The reflection should connect the orbits.")

final_answer = len(orbits) / 2 if is_in_orbit2 and len(orbits) == 2 else 'Error'
if isinstance(final_answer, float): final_answer = int(final_answer)

print(f"\nFinal Equation: {len(orbits)} translation orbits become {final_answer} orbit under full symmetry.")