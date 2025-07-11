import collections

def compose(p1, p2):
    """Computes the composition of two permutations p1 o p2."""
    return [p1[p2[i]] for i in range(len(p1))]

def inverse(p):
    """Computes the inverse of a permutation."""
    inv = [0] * len(p)
    for i, pi in enumerate(p):
        inv[pi] = i
    return inv

def tuple_to_list(p_tuple):
    return list(p_tuple)

def list_to_tuple(p_list):
    return tuple(p_list)

def main():
    n = 3
    
    # All 3! = 6 permutations of {0, 1, 2}
    perms_list = [
        [0, 1, 2], [0, 2, 1], [1, 0, 2],
        [1, 2, 0], [2, 0, 1], [2, 1, 0]
    ]
    perms = {list_to_tuple(p): p for p in perms_list}
    perm_tuples = list(perms.keys())

    # sigma_lht corresponds to the left-hand trefoil
    sigma_lht = [2, 0, 1]

    # Generate all diagrams (X, O) for the left-hand trefoil
    lht_diagrams = set()
    for x_tuple in perm_tuples:
        x = perms[x_tuple]
        o = compose(sigma_lht, x)
        lht_diagrams.add((x_tuple, list_to_tuple(o)))
    
    # Cyclic shift permutations for translations
    shifts = {
        0: [0, 1, 2],
        1: [1, 2, 0],
        2: [2, 0, 1]
    }
    
    # Define symmetry operations
    # Reversal permutation for rotations
    rev = [2, 1, 0]

    def apply_translation(d, r_shift_idx, c_shift_idx):
        x, o = tuple_to_list(d[0]), tuple_to_list(d[1])
        pr_inv = inverse(shifts[r_shift_idx])
        pc = shifts[c_shift_idx]
        
        x_new = compose(pc, compose(x, pr_inv))
        o_new = compose(pc, compose(o, pr_inv))
        return (list_to_tuple(x_new), list_to_tuple(o_new))

    def apply_rotation(d, angle):
        x, o = tuple_to_list(d[0]), tuple_to_list(d[1])
        x_inv, o_inv = inverse(x), inverse(o)
        
        if angle == 0:
            return d
        # Formula for 90-degree rotation on permutations (derived from coordinate transformation)
        if angle == 90:
            x_new = compose(rev, x_inv)
            o_new = compose(rev, o_inv)
        # Formula for 180-degree rotation
        elif angle == 180:
            x_new = compose(rev, compose(x, rev))
            o_new = compose(rev, compose(o, rev))
        # Formula for 270-degree rotation
        elif angle == 270:
            x_new = compose(x_inv, rev)
            o_new = compose(o_inv, rev)
        else:
            return None
        
        return (list_to_tuple(x_new), list_to_tuple(o_new))

    # Check which rotations are valid (map LHT diagrams to LHT diagrams)
    valid_rotations = []
    sample_diagram = next(iter(lht_diagrams))
    for angle in [0, 90, 180, 270]:
        rotated_d = apply_rotation(sample_diagram, angle)
        if rotated_d in lht_diagrams:
            valid_rotations.append(angle)

    # Find the number of orbits (equivalence classes)
    orbits = []
    unclassified = set(lht_diagrams)
    
    while unclassified:
        start_node = unclassified.pop()
        current_orbit = {start_node}
        queue = collections.deque([start_node])
        
        visited = {start_node}
        
        while queue:
            current_d = queue.popleft()
            
            # Apply translations
            for r in range(n):
                for c in range(n):
                    next_d = apply_translation(current_d, r, c)
                    if next_d in lht_diagrams and next_d not in visited:
                        visited.add(next_d)
                        current_orbit.add(next_d)
                        queue.append(next_d)
            
            # Apply valid rotations
            for angle in valid_rotations:
                next_d = apply_rotation(current_d, angle)
                if next_d in lht_diagrams and next_d not in visited:
                    visited.add(next_d)
                    current_orbit.add(next_d)
                    queue.append(next_d)

        orbits.append(current_orbit)
        unclassified.difference_update(current_orbit)

    num_orbits = len(orbits)
    print(f"The {len(lht_diagrams)} minimal grid diagrams for the left-hand trefoil are partitioned into {num_orbits} orbit(s) under translation and rotation.")
    # The actual answer is the number of orbits found
    print(f"The number of unique diagrams is: {num_orbits}")


if __name__ == '__main__':
    main()
