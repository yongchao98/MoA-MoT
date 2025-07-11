import itertools

# Helper functions for permutation operations (0-indexed)
def perm_inverse(p):
    """Computes the inverse of a permutation."""
    inv = [0] * len(p)
    for i, p_i in enumerate(p):
        inv[p_i] = i
    return tuple(inv)

def perm_multiply(p1, p2):
    """Computes the product of two permutations p1 * p2."""
    n = len(p1)
    result = [0] * n
    for i in range(n):
        result[i] = p1[p2[i]]
    return tuple(result)

# Grid parameters
N = 3

# Generate the symmetric group S_3
S3 = list(itertools.permutations(range(N)))

# The 3-cycle corresponding to the left-hand trefoil (0-indexed)
# (1,2,3) -> (3,1,2) becomes (0,1,2) -> (2,0,1)
c2 = (2, 0, 1)

# Generate the 6 grid diagrams for the left-hand trefoil
left_trefoil_diagrams = set()
for sigma in S3:
    tau = perm_multiply(c2, sigma)
    left_trefoil_diagrams.add((sigma, tau))

# Define the D4 group actions on a grid diagram (sigma, tau)
def apply_symmetry(diagram, symmetry_func):
    """Applies a given symmetry function to a diagram."""
    sigma, tau = diagram
    
    # Get the coordinate sets for X's and O's
    x_coords = {(i, sigma[i]) for i in range(N)}
    o_coords = {(i, tau[i]) for i in range(N)}
    
    # Apply the symmetry function to each coordinate
    new_x_coords = {symmetry_func(r, c) for r, c in x_coords}
    new_o_coords = {symmetry_func(r, c) for r, c in o_coords}
    
    # Convert back to (sigma, tau) representation
    new_sigma_list = [0] * N
    new_tau_list = [0] * N
    
    for r, c in new_x_coords:
        new_sigma_list[r] = c
    for r, c in new_o_coords:
        new_tau_list[r] = c
        
    return (tuple(new_sigma_list), tuple(new_tau_list))

# Define the 8 symmetry functions for D4 (0-indexed)
# g(row, col)
d4_symmetries = [
    lambda r, c: (r, c),                             # id
    lambda r, c: (c, N - 1 - r),                     # r90
    lambda r, c: (N - 1 - r, N - 1 - c),             # r180
    lambda r, c: (N - 1 - c, r),                     # r270
    lambda r, c: (r, N - 1 - c),                     # horizontal flip
    lambda r, c: (N - 1 - r, c),                     # vertical flip
    lambda r, c: (c, r),                             # main diagonal flip (y=x)
    lambda r, c: (N - 1 - c, N - 1 - r)              # anti-diagonal flip
]

# Find the number of orbits
num_orbits = 0
remaining_diagrams = set(left_trefoil_diagrams)

# print("The 6 grid diagrams for the left-hand trefoil are:")
# for sigma, tau in sorted(list(remaining_diagrams)):
#     print(f"  (sigma={sigma}, tau={tau})")

while remaining_diagrams:
    num_orbits += 1
    # Pick a diagram to start a new orbit
    start_diagram = remaining_diagrams.pop()
    
    # Find all diagrams in its orbit
    orbit = set()
    for g in d4_symmetries:
        equivalent_diagram = apply_symmetry(start_diagram, g)
        orbit.add(equivalent_diagram)
        
    # Remove the entire orbit from the set of remaining diagrams
    remaining_diagrams.difference_update(orbit)

print("The minimal grid number for the trefoil knot is 3.")
print("There are 6 possible grid diagrams for the left-hand trefoil knot on a 3x3 grid.")
print("When grouping these diagrams by rotational and reflectional symmetry (D4 group),")
print(f"we find that they fall into {num_orbits} distinct classes (orbits).")
print(f"Therefore, there are {num_orbits} grid diagrams for the left-hand trefoil knot up to translation, rotation and reflection.")
