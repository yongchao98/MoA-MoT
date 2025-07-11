import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_problem():
    """
    Finds the proper Euler angle convention equivalent to a given Tait-Bryan rotation.
    """
    # Step 1: Define the initial Tait-Bryan rotation.
    # This is an extrinsic rotation with angles applied in the order X, Y, Z.
    # The scipy convention for this is lowercase 'xyz'.
    tb_angles = [10, 10, 10]
    tb_convention = 'xyz'
    
    # Step 2: Calculate the target rotation matrix.
    target_rotation = Rotation.from_euler(tb_convention, tb_angles, degrees=True)
    target_matrix = target_rotation.as_matrix()

    print("The rotation matrix for the extrinsic Tait-Bryan rotation (X=10°, Y=10°, Z=10°) is:")
    print(np.round(target_matrix, 4))
    print("-" * 30)

    # Step 3: Define the given proper Euler angles and the conventions to test.
    # We assume these are intrinsic rotations (uppercase in scipy), which is common.
    euler_angles = [139.13, 14.11, -141.05]
    proper_conventions = {
        'A': 'XZX',
        'B': 'XYZ', # Not a proper Euler convention, but included for completeness.
        'C': 'YXY',
        'D': 'YZY',
        'E': 'ZYZ',
        'F': 'ZXZ',
    }
    
    found_convention = None
    final_matrix = None

    # Step 4: Iterate through the conventions and compare matrices.
    for key, convention in proper_conventions.items():
        if convention in ['XZX', 'YXY', 'YZY', 'ZYZ', 'ZXZ']: # Only test proper conventions
            test_rotation = Rotation.from_euler(convention, euler_angles, degrees=True)
            test_matrix = test_rotation.as_matrix()
            
            # Check if the generated matrix matches the target matrix
            if np.allclose(target_matrix, test_matrix):
                found_convention = convention
                final_matrix = test_matrix
                break # Exit loop once a match is found

    # Step 5: Display the result.
    if found_convention:
        print(f"A match was found with the '{found_convention}' convention.")
        print("The given Tait-Bryan rotation matrix R is:")
        print("R = R_Z(10.00) R_Y(10.00) R_X(10.00)")
        print(np.round(target_matrix, 4))
        print("\nThis is equivalent to the intrinsic Euler rotation matrix R' for the '{0}' convention:".format(found_convention))
        print("R' = R_{0}({1:.2f}) R_{2}({3:.2f}) R_{4}({5:.2f})".format(
            found_convention[0], euler_angles[0],
            found_convention[1], euler_angles[1],
            found_convention[2], euler_angles[2]
        ))
        print(np.round(final_matrix, 4))
        
        # Determine the letter choice corresponding to the found convention
        for key, value in proper_conventions.items():
            if value == found_convention:
                print(f"\nThis corresponds to answer choice {key}.")
                print(f'<<<{key}>>>')
                break
    else:
        print("No matching convention was found among the given choices.")

solve_rotation_problem()