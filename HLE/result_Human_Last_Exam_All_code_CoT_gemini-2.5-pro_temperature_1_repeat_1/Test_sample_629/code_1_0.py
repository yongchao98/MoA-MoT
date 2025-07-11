import itertools

def perm_inverse(p):
    """Computes the inverse of a permutation."""
    inverse = [0] * len(p)
    for i, val in enumerate(p):
        inverse[val] = i
    return tuple(inverse)

def perm_compose(p1, p2):
    """Computes the composition of two permutations (p1 o p2)."""
    return tuple(p1[p2[i]] for i in range(len(p2)))

def rotate_90(d):
    """Applies a 90-degree clockwise rotation to a grid diagram."""
    X, O = d
    n = len(X)
    w0 = tuple(range(n - 1, -1, -1)) # Reversal permutation, e.g., (2,1,0)
    
    # Based on geometric transformation of coordinates (r,c) -> (c, n-1-r)
    # This leads to X_new = w0 * inv(X) and O_new = w0 * inv(O)
    X_inv = perm_inverse(X)
    O_inv = perm_inverse(O)
    
    new_X = perm_compose(w0, X_inv)
    new_O = perm_compose(w0, O_inv)
    return (new_X, new_O)

def count_lht_grid_diagrams():
    """
    Counts the number of unique minimal grid diagrams for the left-hand trefoil
    knot up to rotation.
    """
    n = 3
    # The characteristic permutation for a left-hand trefoil (1->3, 2->1, 3->2 or 0->2, 1->0, 2->1)
    P_L = (2, 0, 1)

    # 1. Generate all permutations of size n
    perms = list(itertools.permutations(range(n)))

    # 2. Find all grid diagrams (X, O) for the left-hand trefoil
    lht_diagrams = []
    for x_perm in perms:
        for o_perm in perms:
            x_inv = perm_inverse(x_perm)
            norm = perm_compose(x_inv, o_perm)
            if norm == P_L:
                lht_diagrams.append((x_perm, o_perm))

    print(f"Found {len(lht_diagrams)} minimal grid diagrams for the left-hand trefoil knot in total.")
    
    # 3. Group diagrams by rotational equivalence (count orbits)
    unvisited = set(lht_diagrams)
    orbits = []
    while unvisited:
        # Start a new orbit with an arbitrary unvisited diagram
        current_orbit = set()
        q = [unvisited.pop()]
        current_orbit.add(q[0])

        head = 0
        while head < len(q):
            diagram = q[head]
            head += 1
            
            # Find all rotations of the current diagram
            d90 = rotate_90(diagram)
            d180 = rotate_90(d90)
            d270 = rotate_90(d180)
            
            rotations = {d90, d180, d270}
            
            for rot_d in rotations:
                # If a rotation is a valid LHT diagram and not yet in the orbit, add it
                if rot_d in unvisited and rot_d not in current_orbit:
                    current_orbit.add(rot_d)
                    unvisited.remove(rot_d)
                    q.append(rot_d)
        orbits.append(current_orbit)

    print("\nThe diagrams fall into the following orbits under rotation:")
    for i, orbit in enumerate(orbits):
        # To make output readable, convert 0-indexed perms to 1-indexed
        printable_orbit = set()
        for d in orbit:
            p_x = tuple(val + 1 for val in d[0])
            p_o = tuple(val + 1 for val in d[1])
            printable_orbit.add((p_x, p_o))
        print(f"Orbit {i+1}: {printable_orbit}")
    
    num_orbits = len(orbits)
    
    print(f"\nThe equation is the number of distinct orbits found.")
    print(f"Total LHT Diagrams = {len(lht_diagrams)}")
    print(f"Number of diagrams per orbit = {len(lht_diagrams) // num_orbits if num_orbits > 0 else 0}")
    print(f"Number of unique diagrams = {len(lht_diagrams)} / {len(lht_diagrams) // num_orbits if num_orbits > 0 else 0} = {num_orbits}")
    

count_lht_grid_diagrams()