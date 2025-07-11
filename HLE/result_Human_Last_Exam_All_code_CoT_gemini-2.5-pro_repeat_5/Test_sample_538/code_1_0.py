import numpy as np
from itertools import permutations, product

def solve_diffraction_problem():
    """
    Calculates and explains the splitting of Bragg reflections for a rhombohedral
    distortion of a cubic perovskite.
    """
    print("This script determines the number of Bragg reflections for a rhombohedral material (R3m)")
    print("indexed on a pseudocubic cell. The key is to analyze how the peaks from an ideal")
    print("cubic structure split due to the reduction in symmetry.\n")
    print("The rhombohedral distortion occurs along the <111> axis of the pseudocubic cell.")
    print("Planes from a cubic family {hkl} that have different angles with the [111] axis")
    print("will have different d-spacings, leading to separate reflections.\n")
    print("----------------------------------------------------------------------")

    # The unique axis for rhombohedral distortion in a pseudocubic setting
    unique_axis = np.array([1, 1, 1])

    # Define the plane families to be analyzed
    families = {
        "{200}": list(set(permutations([2, 0, 0]))),
        "{220}": list(set(permutations([2, 2, 0]) + permutations([2, -2, 0]))),
        "{222}": list(set(permutations([2, 2, 2]) + permutations([2, 2, -2])))
    }

    total_reflections = 0
    results = {}

    for name, planes in families.items():
        # Store the cosine of the angle for each unique plane orientation
        cos_theta_values = set()

        for plane in planes:
            plane_vec = np.array(plane)
            
            # Calculate the cosine of the angle between the plane normal and the unique axis
            # The absolute value is used as the angle theta and 180-theta are equivalent for splitting
            cos_theta = abs(np.dot(plane_vec, unique_axis) / (np.linalg.norm(plane_vec) * np.linalg.norm(unique_axis)))
            
            # Round to avoid floating point inaccuracies
            cos_theta_values.add(round(cos_theta, 5))

        num_splits = len(cos_theta_values)
        results[name] = num_splits
        total_reflections += num_splits

        print(f"Analysis for the {name} family of planes:")
        print(f"In the rhombohedral structure, the single cubic reflection splits into {num_splits} distinct reflection(s).")
        print("----------------------------------------------------------------------")

    # Final summary and equation
    print("\nSummary of Results:")
    print(f"Number of reflections for {list(results.keys())[0]} = {list(results.values())[0]}")
    print(f"Number of reflections for {list(results.keys())[1]} = {list(results.values())[1]}")
    print(f"Number of reflections for {list(results.keys())[2]} = {list(results.values())[2]}")

    print("\nThe total number of observed Bragg reflections is the sum from each family:")
    
    # Building and printing the final equation string
    equation_parts = [f"{count} ({name})" for name, count in results.items()]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_reflections}")


# Run the main function
solve_diffraction_problem()
<<<5>>>