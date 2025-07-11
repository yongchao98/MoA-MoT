import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_convention():
    """
    Determines the Euler angle convention that produces a rotation equivalent to a given Tait-Bryan rotation.
    """
    # 1. Define the initial Tait-Bryan rotation.
    # The convention is extrinsic X_alpha Y_beta Z_gamma, where alpha=beta=gamma=10 degrees.
    # In scipy, this is represented by the sequence 'XYZ' for extrinsic rotations.
    tb_angles_deg = np.array([10.0, 10.0, 10.0])
    tb_convention = 'XYZ'
    target_rotation = Rotation.from_euler(tb_convention, tb_angles_deg, degrees=True)
    target_matrix = target_rotation.as_matrix()

    print(f"Finding equivalent rotation for extrinsic {tb_convention} with angles:")
    print(f"alpha = {tb_angles_deg[0]}°, beta = {tb_angles_deg[1]}°, gamma = {tb_angles_deg[2]}°\n")

    # 2. Define the Euler angles and conventions to test.
    euler_angles_deg = np.array([139.13, 14.11, -141.05])
    conventions_to_test = {
        'A': 'XZX', 'B': 'XYZ', 'C': 'YXY',
        'D': 'YZY', 'E': 'ZYZ', 'F': 'ZXZ'
    }

    print("Testing against Euler angles:")
    print(f"alpha' = {euler_angles_deg[0]}°, beta' = {euler_angles_deg[1]}°, gamma' = {euler_angles_deg[2]}°\n")
    
    correct_convention_key = None
    correct_convention_str = ""

    # 3. Test each candidate convention (assuming intrinsic rotations, which is standard for Euler angles).
    for key, conv_str in conventions_to_test.items():
        # Scipy uses lowercase for intrinsic conventions (e.g., 'zyz').
        test_rotation = Rotation.from_euler(conv_str.lower(), euler_angles_deg, degrees=True)
        test_matrix = test_rotation.as_matrix()

        # Compare matrices with a tolerance due to rounded input angles.
        is_match = np.allclose(target_matrix, test_matrix, atol=1e-3)
        
        print(f"Testing convention {key} ({conv_str}): {'Match' if is_match else 'No match'}")

        if is_match:
            correct_convention_key = key
            correct_convention_str = conv_str
            
    # 4. Output the final conclusion.
    print("\n---" * 15)
    if correct_convention_key:
        print(f"Conclusion: The correct convention is '{correct_convention_str}' (Choice {correct_convention_key}).\n")
        print("The following two rotations are equivalent:")
        # Outputting the final "equation" with all the numbers
        print(f"1. Extrinsic {tb_convention} rotation with angles (alpha={tb_angles_deg[0]:.2f}°, beta={tb_angles_deg[1]:.2f}°, gamma={tb_angles_deg[2]:.2f}°)")
        print(f"2. Intrinsic {correct_convention_str} rotation with angles (alpha'={euler_angles_deg[0]:.2f}°, beta'={euler_angles_deg[1]:.2f}°, gamma'={euler_angles_deg[2]:.2f}°)")
    else:
        print("No matching convention was found. Please check the input angles.")
    print("---" * 15)

# Execute the function
solve_rotation_convention()