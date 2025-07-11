import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_problem():
    """
    Finds the Euler angle convention that matches a given Tait-Bryan rotation.
    """
    # 1. Define and compute the target rotation.
    # The convention is extrinsic X(alpha)Y(beta)Z(gamma), with alpha=beta=gamma=10 degrees.
    # An extrinsic XYZ rotation is equivalent to an intrinsic ZYX rotation with the same angles.
    # The scipy 'zyx' order means a rotation around z, then the new y, then the newest x.
    # The matrix is R_z(gamma) * R_y(beta) * R_x(alpha).
    # We list angles in the order (gamma, beta, alpha) = (10, 10, 10).
    tait_bryan_angles_deg = [10.0, 10.0, 10.0]
    target_rotation = Rotation.from_euler('zyx', tait_bryan_angles_deg, degrees=True)
    target_matrix = target_rotation.as_matrix()

    # 2. Define the candidate Euler angles and conventions.
    euler_angles_prime_deg = [139.13, 14.11, -141.05]
    
    conventions = {
        'A': 'XZX',
        'B': 'XYZ',
        'C': 'YXY',
        'D': 'YZY',
        'E': 'ZYZ',
        'F': 'ZXZ'
    }

    # 3. Test each convention to find a match.
    correct_choice = None
    for choice, conv_str in conventions.items():
        # For a convention 'ABC', from_euler computes an intrinsic rotation
        # corresponding to R = R_A(angle1) * R_B(angle2) * R_C(angle3).
        candidate_rotation = Rotation.from_euler(conv_str, euler_angles_prime_deg, degrees=True)
        candidate_matrix = candidate_rotation.as_matrix()
        
        # Compare the resulting matrix with the target matrix.
        # A tolerance (atol) is used to account for floating point inaccuracies.
        if np.allclose(target_matrix, candidate_matrix, atol=1e-4):
            correct_choice = choice
            print(f"The equivalent rotation is found using the '{conv_str}' convention (Answer Choice {choice}).")
            print(f"The Euler angles are (alpha', beta', gamma') = ({euler_angles_prime_deg[0]:.2f}°, {euler_angles_prime_deg[1]:.2f}°, {euler_angles_prime_deg[2]:.2f}°).")
            print("The final rotation matrix is:")
            
            # Print each number in the final rotation matrix
            final_matrix = candidate_matrix
            print("[ {:.8f}  {:.8f}  {:.8f} ]".format(final_matrix[0,0], final_matrix[0,1], final_matrix[0,2]))
            print("[ {:.8f}  {:.8f}  {:.8f} ]".format(final_matrix[1,0], final_matrix[1,1], final_matrix[1,2]))
            print("[ {:.8f}  {:.8f}  {:.8f} ]".format(final_matrix[2,0], final_matrix[2,1], final_matrix[2,2]))
            break
            
    if correct_choice is None:
        print("No matching convention was found in the given options.")

solve_rotation_problem()