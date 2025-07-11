import itertools

def perm_to_dict(p_tuple, n=3):
    """Converts a tuple representation of a permutation into a dictionary."""
    d = {i: i for i in range(1, n + 1)}
    for cycle in p_tuple:
        if len(cycle) > 1:
            for i in range(len(cycle)):
                d[cycle[i]] = cycle[(i + 1) % len(cycle)]
    return d

def dict_to_tuple(p_dict, n=3):
    """Converts a dict representation of a permutation back to a tuple."""
    tuples = []
    processed = [False] * (n + 1)
    for i in range(1, n + 1):
        if processed[i]:
            continue
        cycle = []
        j = i
        while not processed[j]:
            processed[j] = True
            cycle.append(j)
            j = p_dict.get(j, j)
            if j == i:
                break
        if len(cycle) > 1:
            tuples.append(tuple(cycle))
    # Sort for canonical representation
    return tuple(sorted(tuples))

def multiply_perms(p2, p1):
    """Multiplies two permutations p2 * p1. Input are dicts."""
    res_dict = {}
    for i in p1:
        res_dict[i] = p2.get(p1.get(i, i), p1.get(i, i))
    return res_dict

def invert_perm(p):
    """Inverts a permutation. Input is a dict."""
    return {v: k for k, v in p.items()}

# Define S3 permutations as tuples
s3_tuples = [
    (),  # identity
    ((1, 2),),
    ((1, 3),),
    ((2, 3),),
    ((1, 2, 3),),
    ((1, 3, 2),)
]

# Convert S3 to dictionary representation for easier computation
s3_dicts = [perm_to_dict(p) for p in s3_tuples]

# pi_L for Left-Hand Trefoil
pi_L_tuple = ((1, 3, 2),)
pi_L_dict = perm_to_dict(pi_L_tuple)

# Generate the 6 LHT diagrams (X, O) where O * X^-1 = pi_L, so O = pi_L * X
lht_diagrams = []
print("Generating 6 diagrams for the Left-Hand Trefoil:")
for x_dict in s3_dicts:
    o_dict = multiply_perms(pi_L_dict, x_dict)
    # Convert back to tuple for display and storage
    x_tuple = dict_to_tuple(x_dict)
    o_tuple = dict_to_tuple(o_dict)
    lht_diagrams.append((x_tuple, o_tuple))
    print(f"X = {x_tuple}, O = {o_tuple}")

# Define the translation operator
c_tuple = ((1, 2, 3),)
c_dict = perm_to_dict(c_tuple)
c_inv_dict = invert_perm(c_dict)

def apply_translation(diagram_tuple):
    x_tuple, o_tuple = diagram_tuple
    x_dict, o_dict = perm_to_dict(x_tuple), perm_to_dict(o_tuple)
    
    # new_X = c * X * c^-1
    new_x_dict = multiply_perms(c_dict, multiply_perms(x_dict, c_inv_dict))
    
    # new_O = c * O * c^-1
    new_o_dict = multiply_perms(c_dict, multiply_perms(o_dict, c_inv_dict))
    
    return (dict_to_tuple(new_x_dict), dict_to_tuple(new_o_dict))

# Find the orbits
print("\nCalculating orbits under translation:")
remaining_diagrams = set(lht_diagrams)
orbit_count = 0
orbits = []

while remaining_diagrams:
    orbit_count += 1
    orbit = set()
    # Pick a starting diagram
    start_diagram = remaining_diagrams.pop()
    orbit.add(start_diagram)
    
    # Generate the rest of the orbit
    current_diagram = start_diagram
    while True:
        next_diagram = apply_translation(current_diagram)
        if next_diagram in orbit:
            break
        orbit.add(next_diagram)
        current_diagram = next_diagram
        
    orbits.append(orbit)
    remaining_diagrams -= orbit

print(f"\nFound {orbit_count} orbits.")
for i, orbit in enumerate(orbits):
    print(f"Orbit {i+1}: {orbit}")
    
print(f"\nFinal Answer: The number of distinct grid diagrams is the number of orbits.")
print(f"1 + 1 + 1 + 1 = 4")
