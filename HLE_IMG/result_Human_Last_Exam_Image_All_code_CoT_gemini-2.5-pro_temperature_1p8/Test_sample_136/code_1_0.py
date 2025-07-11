import itertools

def generate_colorings():
    """Generates all colorings of 6 beads with 1 Blue and 1 Green."""
    beads = list(range(6))
    # In C(n,k), n=6 beads, k=2 colored beads. C(6,2) = 15.
    # But Blue/Green are distinct, so it's P(6,2) = 6*5 = 30.
    colorings = set()
    for b_pos in beads:
        remaining_beads = [p for p in beads if p != b_pos]
        for g_pos in remaining_beads:
            c = ['P'] * 6
            c[b_pos] = 'B'
            c[g_pos] = 'G'
            colorings.add(tuple(c))
    return colorings

def apply_permutation(coloring, perm):
    """
    Applies a permutation to a coloring.
    The new coloring at position i gets its color from the old coloring
    at position p_inv[i], where p_inv is the inverse of the permutation.
    """
    perm_len = len(perm)
    p_inv = [0] * perm_len
    for i, p_i in enumerate(perm):
        p_inv[p_i] = i
    
    new_coloring = [None] * perm_len
    for i in range(perm_len):
        new_coloring[i] = coloring[p_inv[i]]
    return tuple(new_coloring)

def count_orbits_burnside(colorings, group):
    """Counts the number of distinct orbits using Burnside's Lemma."""
    total_fixed_points = 0
    for g in group:
        fixed_count = 0
        for c in colorings:
            if apply_permutation(c, g) == c:
                fixed_count += 1
        total_fixed_points += fixed_count
    
    num_orbits = total_fixed_points / len(group)
    return int(num_orbits)

def solve():
    """Identifies the symmetry group by calculating orbit counts."""
    print("Plan: Use Burnside's Lemma to find the symmetry group that produces 5 equivalence classes.")
    
    # 1. Setup colorings and group generators
    n_beads = 6
    colorings = generate_colorings()
    identity = tuple(range(n_beads))
    
    # Generator for rotation (C6)
    r = tuple((i + 1) % n_beads for i in range(n_beads)) # (1, 2, 3, 4, 5, 0)
    
    # Generator for reflection (to generate D6)
    s = tuple( (n_beads - i) % n_beads for i in range(n_beads)) # (0, 5, 4, 3, 2, 1)

    # 2. Generate C6 group (rotations)
    c6_elements = {identity}
    current_perm = r
    for _ in range(n_beads - 1):
        c6_elements.add(current_perm)
        # Compose r with current_perm: new_p(i) = r(current_p(i))
        current_perm = tuple(r[i] for i in current_perm)
        
    # 3. Generate D6 group (rotations and reflections)
    reflections = set()
    for p_r in c6_elements:
        # Compose p_r with s: new_p(i) = p_r(s(i))
        new_reflection = tuple(p_r[s[i]] for i in range(n_beads))
        reflections.add(new_reflection)
    d6_elements = c6_elements.union(reflections)

    # 4. Calculate number of orbits for each group
    num_c6_orbits = count_orbits_burnside(colorings, c6_elements)
    num_d6_orbits = count_orbits_burnside(colorings, d6_elements)
    
    # 5. Print results and conclusion
    print("\n--- Calculation Results ---")
    print(f"Total number of unique colorings (1 Blue, 1 Green): {len(colorings)}")
    print(f"Number of orbits produced by the C6 (rotations only) group: {num_c6_orbits}")
    print(f"Number of orbits produced by the D6 (rotations & reflections) group: {num_d6_orbits}")
    
    print("\n--- Conclusion ---")
    print("The problem image shows 5 distinct rows, where each row is an equivalence class.")
    print(f"Our calculation shows that the C6 group produces {num_c6_orbits} classes, which matches the image.")
    print("The D6 group produces a different number of classes.")
    print("Therefore, the group of symmetries is the cyclic group C6.")
    
    print("\nThe minimal generator for C6 on a hexagon is the smallest positive rotation.")
    print("Minimal generator: rotation by 60 degrees")

solve()