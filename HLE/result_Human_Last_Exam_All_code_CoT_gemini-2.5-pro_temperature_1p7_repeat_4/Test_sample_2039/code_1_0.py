import numpy as np
from scipy.spatial.transform import Rotation

def find_equivalent_rotation_convention():
    """
    Calculates the rotation matrix for a given Tait-Bryan rotation and finds which
    proper Euler angle convention with given angles produces the same rotation.
    """
    # Define the extrinsic Tait-Bryan angles in degrees
    alpha = 10.0
    beta = 10.0
    gamma = 10.0
    tb_angles = [alpha, beta, gamma]
    
    # Define the equivalent proper Euler angles in degrees
    alpha_p = 139.13
    beta_p = 14.11
    gamma_p = -141.05
    pe_angles = [alpha_p, beta_p, gamma_p]
    
    # The initial rotation is given by an extrinsic X-Y-Z convention.
    # The scipy 'XYZ' convention calculates the rotation matrix as R_Z(gamma) * R_Y(beta) * R_X(alpha).
    target_rotation = Rotation.from_euler('XYZ', tb_angles, degrees=True)
    target_matrix = target_rotation.as_matrix()

    # Define the proper Euler conventions from the answer choices.
    # We assume these are intrinsic conventions (lowercase in scipy).
    candidate_conventions = {
        'A': 'xzx',
        'C': 'yxy',
        'D': 'yzy',
        'E': 'zyz',
        'F': 'zxz'
    }

    print("Comparing the rotation from extrinsic Tait-Bryan angles:")
    print(f"alpha = {alpha} deg, beta = {beta} deg, gamma = {gamma} deg (around X, Y, Z axes)")
    print("\nwith rotations from various proper Euler angle conventions with angles:")
    print(f"alpha' = {alpha_p} deg, beta' = {beta_p} deg, gamma' = {gamma_p} deg\n")

    # Iterate through the candidate conventions and find the match
    for key, conv_str in candidate_conventions.items():
        candidate_rotation = Rotation.from_euler(conv_str, pe_angles, degrees=True)
        candidate_matrix = candidate_rotation.as_matrix()
        
        # Check if the matrices are close enough to be considered identical
        if np.allclose(target_matrix, candidate_matrix, atol=1e-4):
            print("--------------------------------------------------")
            print(f"SUCCESS: A match was found!")
            print(f"The rotation from extrinsic Tait-Bryan convention 'XYZ' with angles ({alpha}, {beta}, {gamma}) is equivalent to")
            print(f"the rotation from intrinsic proper Euler convention '{conv_str.upper()}' with angles ({alpha_p}, {beta_p}, {gamma_p}).")
            print(f"\nThis corresponds to answer choice {key}.")
            print("--------------------------------------------------")
            return f"<<<{key}>>>"
    
    print("Could not find a matching convention among the choices.")
    return None

# Execute the function and print the final result
find_equivalent_rotation_convention()
