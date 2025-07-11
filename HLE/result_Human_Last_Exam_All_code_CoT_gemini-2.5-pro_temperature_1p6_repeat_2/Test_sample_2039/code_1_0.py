import numpy as np
from scipy.spatial.transform import Rotation

def find_rotation_convention():
    """
    Finds the correct Euler angle convention that matches a given rotation.
    """
    # Step 1: Define the initial extrinsic rotation and compute its matrix.
    # The convention is extrinsic XYZ with angles alpha=10, beta=10, gamma=10 degrees.
    initial_angles = [10.0, 10.0, 10.0]
    # In scipy, 'XYZ' denotes an extrinsic rotation. The matrix is R_z(10) * R_y(10) * R_x(10).
    r_initial = Rotation.from_euler('XYZ', initial_angles, degrees=True)
    R_initial_matrix = r_initial.as_matrix()

    print(f"The reference rotation matrix from extrinsic XYZ with angles ({initial_angles[0]}, {initial_angles[1]}, {initial_angles[2]}) is:")
    print(R_initial_matrix)
    print("-" * 50)

    # Step 2: Define the target Euler angles to test.
    target_angles = [139.13, 14.11, -141.05]
    print(f"Testing conventions with Euler angles ({target_angles[0]}, {target_angles[1]}, {target_angles[2]})...\n")


    # Step 3: Iterate through each candidate convention and compare matrices.
    # The term "Euler angles" usually implies an intrinsic rotation sequence,
    # which is represented by lowercase letters in scipy.
    answer_choices = {
        "A": "XZX",
        "B": "XYZ",
        "C": "YXY",
        "D": "YZY",
        "E": "ZYZ",
        "F": "ZXZ",
    }

    correct_convention = None
    for key, conv_str in answer_choices.items():
        # Test the intrinsic convention (e.g., 'zyz')
        r_test = Rotation.from_euler(conv_str.lower(), target_angles, degrees=True)
        R_test_matrix = r_test.as_matrix()

        # Compare the test matrix with the reference matrix
        is_match = np.allclose(R_initial_matrix, R_test_matrix, atol=1e-4)

        print(f"Testing convention {key} ({conv_str}):")
        if is_match:
            print(f"MATCH FOUND for convention {conv_str}.")
            print("The matrix for this convention is:")
            print(R_test_matrix)
            correct_convention = key
            break  # Stop once the match is found
        else:
            print("No match.")
        print()

    if not correct_convention:
        print("No matching convention was found among the given options.")

find_rotation_convention()