import numpy as np

def calculate_peak_splitting(h, k, l, distortion_axis):
    """
    Calculates the number of split Bragg reflections for a {hkl} family.
    The splitting is due to a structural distortion from cubic to a lower symmetry
    phase (e.g., rhombohedral) along a specific axis.

    Args:
        h, k, l (int): The Miller indices of the plane family in pseudocubic notation.
        distortion_axis (np.array): The vector representing the unique distortion axis.

    Returns:
        int: The number of observable Bragg reflections.
    """
    # Generate a set of representative vectors for the cubic {hkl} family.
    # For example, for {200}, the representative directions are [2,0,0], [0,2,0], [0,0,2].
    family_vectors = set()
    # Use permutations to get planes like (h,k,l), (k,h,l), etc.
    for p in set(np.abs(p) for p in np.unique(np.array(list(map(list,set(p for p in __import__('itertools').permutations([h,k,l]))))),axis=0)):
        # Use product to get sign variations like (h,k,l), (h,-k,l), etc.
        # We only need to check one from each pair of opposite vectors, e.g., (h,k,l) and (-h,-k,-l)
        for s in __import__('itertools').product([-1, 1], repeat=3):
            vec = tuple(sorted(np.abs(p * s)))
            family_vectors.add(vec)

    # Convert to a list of numpy arrays for vector math
    family_vectors = [np.array(v) for v in family_vectors if np.any(v)]

    d_norm_sq = np.dot(distortion_axis, distortion_axis)
    
    # Use a set to store the unique values of cos^2(angle)
    # The number of unique values is the number of split peaks.
    cos_sq_values = set()
    
    for v in family_vectors:
        v_norm_sq = np.dot(v, v)
        if v_norm_sq == 0:
            continue
            
        # Calculate the square of the cosine of the angle between the plane normal and the distortion axis
        dot_product = np.dot(v, distortion_axis)
        cos_sq = (dot_product**2) / (v_norm_sq * d_norm_sq)
        
        # Round to a few decimal places to avoid floating-point inaccuracies
        cos_sq_values.add(round(cos_sq, 5))
        
    return len(cos_sq_values)

def main():
    """
    Main function to solve the problem and print the results.
    """
    # The distortion from cubic to rhombohedral (R3m) occurs along the <111> direction.
    distortion_axis = np.array([1, 1, 1])

    # Families of planes to analyze (in pseudocubic indexing)
    plane_families = {
        "{200}": (2, 0, 0),
        "{220}": (2, 2, 0),
        "{222}": (2, 2, 2)
    }

    print("Calculating the number of Bragg reflections for a rhombohedral (R3m) material...\n")

    results = {}
    for name, indices in plane_families.items():
        h, k, l = indices
        num_peaks = calculate_peak_splitting(h, k, l, distortion_axis)
        results[name] = num_peaks

    # Output the results
    n1 = results["{200}"]
    n2 = results["{220}"]
    n3 = results["{222}"]
    total_peaks = n1 + n2 + n3

    print(f"For the {{200}} family of planes, the number of reflections is: {n1}")
    print(f"For the {{220}} family of planes, the number of reflections is: {n2}")
    print(f"For the {{222}} family of planes, the number of reflections is: {n3}\n")
    print("The final equation for the total number of reflections from these families is:")
    print(f"Total Reflections = {n1} + {n2} + {n3} = {total_peaks}")

if __name__ == "__main__":
    main()
<<<5>>>