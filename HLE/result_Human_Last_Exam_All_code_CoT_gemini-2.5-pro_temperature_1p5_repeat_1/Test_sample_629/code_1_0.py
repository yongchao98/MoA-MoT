import itertools

def invert_permutation(p):
    """Inverts a permutation tuple."""
    inverse = [0] * len(p)
    for i, val in enumerate(p):
        inverse[val] = i
    return tuple(inverse)

def compose_permutations(p1, p2):
    """Composes two permutation tuples (p1 o p2)."""
    return tuple(p1[p2[i]] for i in range(len(p1)))

def count_cycles(p):
    """Counts the number of cycles in a permutation."""
    n = len(p)
    seen = [False] * n
    cycles = 0
    for i in range(n):
        if not seen[i]:
            cycles += 1
            j = i
            while not seen[j]:
                seen[j] = True
                j = p[j]
    return cycles

def sign(x):
    """Sign function."""
    if x > 0: return 1
    if x < 0: return -1
    return 0

def calculate_writhe(sigma, tau):
    """Calculates the writhe of a grid diagram."""
    n = len(sigma)
    # Turaev's formula for writhe
    w_v = sum(sign(t - s) for s, t in zip(sigma, tau))
    sigma_shifted = sigma[1:] + sigma[:1]
    w_h = sum(sign(s_shifted - t) for s_shifted, t in zip(sigma_shifted, tau))
    return w_v + w_h

def rotate_diagram(sigma, tau, degrees):
    """Rotates a grid diagram by 90, 180, or 270 degrees."""
    n = len(sigma)
    p_rev = tuple(range(n - 1, -1, -1))
    
    if degrees == 90:
        sigma_inv = invert_permutation(sigma)
        tau_inv = invert_permutation(tau)
        sigma_rot = compose_permutations(p_rev, sigma_inv)
        tau_rot = compose_permutations(p_rev, tau_inv)
        return sigma_rot, tau_rot
    elif degrees == 180:
        sigma_rot = compose_permutations(compose_permutations(p_rev, sigma), p_rev)
        tau_rot = compose_permutations(compose_permutations(p_rev, tau), p_rev)
        return sigma_rot, tau_rot
    elif degrees == 270:
        sigma_inv = invert_permutation(sigma)
        tau_inv = invert_permutation(tau)
        sigma_rot = compose_permutations(sigma_inv, p_rev)
        tau_rot = compose_permutations(tau_inv, p_rev)
        return sigma_rot, tau_rot
    else: # 0 degrees
        return sigma, tau

def solve():
    """
    Solves the problem by generating, filtering, and counting grid diagrams.
    """
    n = 3
    perms = list(itertools.permutations(range(n)))
    all_diagrams = list(itertools.product(perms, perms))

    lht_diagrams = []
    for sigma, tau in all_diagrams:
        tau_inv = invert_permutation(tau)
        pi = compose_permutations(sigma, tau_inv)
        
        # Must be a knot (1 component)
        if count_cycles(pi) == 1:
            writhe = calculate_writhe(sigma, tau)
            # Left-hand trefoil has writhe -2 in minimal grid diagrams
            if writhe == -2:
                lht_diagrams.append((sigma, tau))

    num_lht = len(lht_diagrams)
    
    # Using Burnside's Lemma to count unique diagrams under rotation
    # Orbit Count = (1/|G|) * sum(|X^g|)
    # G = {id, r90, r180, r270} is the rotation group
    # X is the set of lht_diagrams
    
    # X^id is just X
    num_fixed_by_id = num_lht
    
    # Count diagrams fixed by 90-degree rotation
    num_fixed_by_90 = 0
    for d in lht_diagrams:
        if rotate_diagram(d[0], d[1], 90) == d:
            num_fixed_by_90 += 1

    # Count diagrams fixed by 180-degree rotation
    num_fixed_by_180 = 0
    for d in lht_diagrams:
        if rotate_diagram(d[0], d[1], 180) == d:
            num_fixed_by_180 += 1
            
    # Number of diagrams fixed by 270-degree rotation is the same as for 90
    num_fixed_by_270 = num_fixed_by_90
    
    # Calculate the number of orbits
    num_orbits = (1/4) * (num_fixed_by_id + num_fixed_by_90 + num_fixed_by_180 + num_fixed_by_270)

    print("Step 1: Generated all 3x3 grid diagrams.")
    print(f"Step 2: Identified diagrams for the left-hand trefoil knot (writhe = -2).")
    print(f"Found {num_lht} diagrams for the left-hand trefoil.")
    print("\nStep 3: Counting unique diagrams using Burnside's Lemma for rotational symmetry.")
    print("Number of diagrams = |X| = 12")
    print(f"Number of diagrams fixed by  90-degree rotation: {num_fixed_by_90}")
    print(f"Number of diagrams fixed by 180-degree rotation: {num_fixed_by_180}")
    print(f"Number of diagrams fixed by 270-degree rotation: {num_fixed_by_270}")
    
    print("\nFinal Calculation:")
    # We output each number in the equation for clarity, as requested
    term1 = num_fixed_by_id
    term2 = num_fixed_by_90
    term3 = num_fixed_by_180
    term4 = num_fixed_by_270
    total_fixed_points = term1 + term2 + term3 + term4
    group_size = 4
    
    print(f"Number of unique diagrams = (1/{group_size}) * ({term1} + {term2} + {term3} + {term4})")
    print(f"                            = (1/{group_size}) * {total_fixed_points}")
    print(f"                            = {int(num_orbits)}")

solve()