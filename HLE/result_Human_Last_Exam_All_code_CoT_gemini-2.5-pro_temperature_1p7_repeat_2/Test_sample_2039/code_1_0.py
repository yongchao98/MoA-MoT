import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_problem():
    """
    Determines the correct proper Euler angle convention for a given rotation.
    """
    # Define the initial extrinsic Tait-Bryan rotation
    alpha, beta, gamma = 10.0, 10.0, 10.0
    tait_bryan_angles = [alpha, beta, gamma]
    
    # Calculate the target rotation matrix.
    # The 'xyz' convention in scipy.spatial.transform.Rotation.from_euler
    # corresponds to extrinsic rotations around X, then Y, then Z.
    target_rotation = Rotation.from_euler('xyz', tait_bryan_angles, degrees=True)
    target_matrix = target_rotation.as_matrix()

    print(f"The task is to find the proper Euler angle convention that is equivalent to an extrinsic XYZ rotation.")
    print(f"The equation for the initial extrinsic rotation is R = Rz({gamma:.2f}°) Ry({beta:.2f}°) Rx({alpha:.2f}°)")
    print("\nResulting target rotation matrix:")
    print(np.round(target_matrix, 4))
    
    # Define the equivalent Euler angles and the conventions to test
    equivalent_angles = [139.13, 14.11, -141.05]
    alpha_p, beta_p, gamma_p = equivalent_angles
    
    # Note: Proper Euler conventions have the first and third axes the same.
    # The convention strings in scipy ('xzx', 'yxy', etc.) represent intrinsic rotations.
    conventions_to_test = {
        'A': 'xzx',
        'B': 'xyz', # This is Tait-Bryan, not proper Euler, but we can check.
        'C': 'yxy',
        'D': 'yzy',
        'E': 'zyz',
        'F': 'zxz'
    }

    print("\n--- Testing Potential Conventions ---")

    matching_convention = None
    
    for choice, conv_str in conventions_to_test.items():
        test_rotation = Rotation.from_euler(conv_str, equivalent_angles, degrees=True)
        test_matrix = test_rotation.as_matrix()
        
        # Check if the calculated matrix is close to the target matrix
        if np.allclose(target_matrix, test_matrix):
            matching_convention = conv_str.upper()
            
            print(f"\nFound a match with the {matching_convention} convention.")
            print(f"The equation for the equivalent intrinsic rotation is R = R_{conv_str[0].upper()}({alpha_p:.2f}°) R_{conv_str[1].upper()}({beta_p:.2f}°) R_{conv_str[2].upper()}({gamma_p:.2f}°)")
            print("\nResulting rotation matrix from this convention:")
            print(np.round(test_matrix, 4))
            print("\nThe two matrices are identical, so this is the correct convention.")
            final_answer = f"<<<{choice}>>>"
            break
            
    if not matching_convention:
        print("\nCould not find a matching convention among the choices.")
        final_answer = "<<<No match found>>>"
        
    return final_answer

# Run the solver and print the final answer
final_answer_result = solve_rotation_problem()
print(final_answer_result)