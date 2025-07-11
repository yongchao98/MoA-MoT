import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_convention():
    """
    Finds the correct Euler angle convention that matches a given Tait-Bryan rotation.
    """
    # Define the initial extrinsic Tait-Bryan rotation angles in degrees.
    # The convention is X-Y-Z, with rotations applied in that order around fixed axes.
    tb_angles = [10.0, 10.0, 10.0]
    tb_convention = 'xyz' # 'xyz' (lowercase) for extrinsic rotations in scipy

    # Calculate the rotation object for the initial Tait-Bryan rotation.
    target_rotation = Rotation.from_euler(tb_convention, tb_angles, degrees=True)

    # Define the given proper Euler angles in degrees.
    euler_angles = [139.13, 14.11, -141.05]

    # Define the possible Euler angle conventions from the answer choices.
    # Proper Euler conventions are intrinsic and denoted by uppercase letters in scipy.
    conventions_to_test = {
        "A": "XZX",
        "B": "XYZ", # Tait-Bryan, not proper Euler
        "C": "YXY",
        "D": "YZY",
        "E": "ZYZ",
        "F": "ZXZ"
    }

    found_convention_key = None
    found_convention_name = ""

    # Iterate through each convention and check if it matches the target rotation.
    for key, name in conventions_to_test.items():
        try:
            # Calculate the rotation object for the current convention and angles.
            test_rotation = Rotation.from_euler(name, euler_angles, degrees=True)
            
            # Compare the two rotations. np.allclose is used to check if the
            # underlying rotation matrices are nearly identical.
            if np.allclose(target_rotation.as_matrix(), test_rotation.as_matrix(), atol=1e-4):
                found_convention_key = key
                found_convention_name = name
                break
        except ValueError:
            # Handle any potential errors with invalid convention strings
            print(f"Could not process convention: {name}")

    # Output the result
    if found_convention_key:
        print(f"The initial Tait-Bryan rotation using the {tb_convention.upper()} convention is:")
        print(f"alpha = {tb_angles[0]}°, beta = {tb_angles[1]}°, gamma = {tb_angles[2]}°")
        print("\nThis is equivalent to the proper Euler rotation using the following convention and angles:")
        print(f"Convention: {found_convention_name} (Choice {found_convention_key})")
        print(f"alpha' = {euler_angles[0]}°, beta' = {euler_angles[1]}°, gamma' = {euler_angles[2]}°")
        print(f"\n<<<E>>>")
    else:
        print("No matching convention was found among the choices.")

solve_rotation_convention()