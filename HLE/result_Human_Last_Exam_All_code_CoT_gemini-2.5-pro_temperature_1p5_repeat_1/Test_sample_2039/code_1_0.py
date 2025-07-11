import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_problem():
    """
    Finds the equivalent proper Euler angle convention for a given Tait-Bryan rotation.
    """
    # Step 1: Define the initial extrinsic rotation (Tait-Bryan angles).
    # Convention: XYZ (extrinsic), Angles: a=10, b=10, g=10 degrees.
    initial_convention = 'XYZ'
    initial_angles_deg = [10, 10, 10]

    # Step 2: Create the reference rotation object and its matrix representation.
    # The 'XYZ' string with uppercase letters specifies extrinsic rotations.
    r_ref = Rotation.from_euler(initial_convention, initial_angles_deg, degrees=True)
    R_ref = r_ref.as_matrix()

    # The equivalent proper Euler angles to be tested.
    euler_angles_prime_deg = [139.13, 14.11, -141.05]

    # Step 3: Define the candidate proper Euler conventions from the answer choices.
    candidate_conventions = {
        'A': 'XZX',
        'B': 'XYZ', # This is a Tait-Bryan convention, not proper Euler, but included for completeness.
        'C': 'YXY',
        'D': 'YZY',
        'E': 'ZYZ',
        'F': 'ZXZ'
    }

    print(f"Finding the proper Euler angle convention equivalent to the extrinsic rotation:")
    print(f"R({initial_convention}, α={initial_angles_deg[0]}°, β={initial_angles_deg[1]}°, γ={initial_angles_deg[2]}°)\n")

    matching_convention_key = None
    
    # Step 4 & 5: Iterate through candidates, create rotation matrices, and compare.
    for key, conv_str in candidate_conventions.items():
        # Proper Euler conventions are 'ZXZ', 'XYX', 'YZY', 'ZYZ', 'XZX', 'YXY'.
        # We only test the ones provided in the multiple-choice question.
        if key in ['A', 'C', 'D', 'E', 'F']:
            r_test = Rotation.from_euler(conv_str, euler_angles_prime_deg, degrees=True)
            R_test = r_test.as_matrix()

            # np.allclose checks if two arrays are element-wise equal within a tolerance.
            if np.allclose(R_ref, R_test, atol=1e-4):
                matching_convention_key = key
                # Step 6: A match is found, print the detailed result.
                print(f"--- MATCH FOUND ---")
                print(f"The rotation from extrinsic convention '{initial_convention}' with angles α={initial_angles_deg[0]}°, β={initial_angles_deg[1]}°, γ={initial_angles_deg[2]}°")
                print("is equivalent to")
                print(f"the rotation from proper Euler convention '{conv_str}' ({key}) with angles α'={euler_angles_prime_deg[0]}°, β'={euler_angles_prime_deg[1]}°, γ'={euler_angles_prime_deg[2]}°.")
                break # Exit the loop as we found our answer.

    if not matching_convention_key:
        print("\nNo matching convention was found among the given options.")

# Execute the function to solve the problem
solve_rotation_problem()
<<<F>>>