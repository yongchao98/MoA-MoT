import itertools

def perm_compose(p1, p2):
    """Computes the composition of two permutations."""
    return tuple(p1[i] for i in p2)

def perm_inverse(p):
    """Computes the inverse of a permutation."""
    n = len(p)
    inv = [0] * n
    for i, val in enumerate(p):
        inv[val] = i
    return tuple(inv)

def apply_row_translation(diagram, k, n):
    """Applies row translation to a diagram."""
    x, o = diagram
    return (tuple((val + k) % n for val in x), tuple((val + k) % n for val in o))

def apply_col_translation(diagram, k, n):
    """Applies column translation to a diagram."""
    x, o = diagram
    return (tuple(x[(i - k + n) % n] for i in range(n)), tuple(o[(i - k + n) % n] for i in range(n)))

def apply_180_rotation(diagram, n):
    """Applies 180-degree rotation to a diagram."""
    x, o = diagram
    rho = tuple(range(n - 1, -1, -1))
    x_new = perm_compose(rho, perm_compose(x, rho))
    o_new = perm_compose(rho, perm_compose(o, rho))
    return (x_new, o_new)

def main():
    """
    Counts the number of unique minimal grid diagrams for the left-hand trefoil
    up to translation and 180-degree rotation.
    """
    n = 3
    # The permutation cycle corresponding to the left-hand trefoil
    c_lht = (2, 0, 1)

    # Generate all 6 minimal grid diagrams for the left-hand trefoil
    all_perms = list(itertools.permutations(range(n)))
    lht_diagrams = set()
    for o_perm in all_perms:
        # Knot type is determined by O_inv * X. So X = O * C.
        x_perm = perm_compose(o_perm, c_lht)
        lht_diagrams.add((x_perm, o_perm))

    # Group diagrams into equivalence classes (orbits)
    orbits = []
    unclassified_diagrams = set(lht_diagrams)

    while unclassified_diagrams:
        # Start a new orbit with a diagram
        seed = unclassified_diagrams.pop()
        current_orbit = {seed}
        queue = [seed]

        # Explore all diagrams reachable from the seed
        while queue:
            current_diagram = queue.pop(0)

            # Apply all allowed transformations
            for k in range(1, n):
                # Row and column translations
                g_tr = apply_row_translation(current_diagram, k, n)
                g_tc = apply_col_translation(current_diagram, k, n)
                
                if g_tr not in current_orbit:
                    current_orbit.add(g_tr)
                    queue.append(g_tr)
                if g_tc not in current_orbit:
                    current_orbit.add(g_tc)
                    queue.append(g_tc)

            # 180-degree rotation
            g_r180 = apply_180_rotation(current_diagram, n)
            if g_r180 not in current_orbit:
                current_orbit.add(g_r180)
                queue.append(g_r180)

        # The orbit within the set of LHT diagrams is the intersection
        orbit_in_lht_set = current_orbit.intersection(lht_diagrams)
        orbits.append(orbit_in_lht_set)
        unclassified_diagrams -= orbit_in_lht_set

    print("The number of unique minimal grid diagrams for the left-hand trefoil knot is:")
    print(len(orbits))

if __name__ == "__main__":
    main()