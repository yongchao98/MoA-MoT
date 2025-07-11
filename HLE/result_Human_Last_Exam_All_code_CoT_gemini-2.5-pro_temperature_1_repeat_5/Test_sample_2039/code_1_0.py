import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_convention():
    """
    Finds the correct Euler angle convention that matches a given Tait-Bryan rotation.
    """
    # 1. Define the initial Tait-Bryan rotation and compute the target matrix.
    # Convention is X-Y-Z extrinsic, with alpha=beta=gamma=10 degrees.
    # In scipy, this is represented by the string 'xyz'.
    tb_angles = [10, 10, 10]
    r_target = Rotation.from_euler('xyz', tb_angles, degrees=True)
    m_target = r_target.as_matrix()

    # 2. Define the Euler angles to test.
    euler_angles = [139.13, 14.11, -141.05]

    # 3. Define the possible proper Euler conventions from the answer choices.
    conventions = {
        "A": "xzx",
        "B": "xyz", # Tait-Bryan, not proper Euler, but included for completeness.
        "C": "yxy",
        "D": "yzy",
        "E": "zyz",
        "F": "zxz"
    }

    # 4. Iterate through the conventions and check for a match.
    found_convention_key = None
    for key, convention_str in conventions.items():
        r_test = Rotation.from_euler(convention_str, euler_angles, degrees=True)
        m_test = r_test.as_matrix()
        
        # Compare the test matrix with the target matrix.
        if np.allclose(m_target, m_test, atol=1e-4):
            found_convention_key = key
            break

    # 5. Print the result.
    if found_convention_key:
        found_convention_str = conventions[found_convention_key].upper()
        print(f"The target rotation from the Tait-Bryan angles ({tb_angles[0]}°, {tb_angles[1]}°, {tb_angles[2]}°) using the XYZ convention is equivalent to the rotation from:")
        print(f"The Euler angles ({euler_angles[0]}°, {euler_angles[1]}°, {euler_angles[2]}°) using the {found_convention_str} convention.")
        print(f"\nThe correct answer choice is {found_convention_key}.")
    else:
        print("No matching convention was found among the options.")

    # Return the final answer key.
    return found_convention_key

# Run the solver and print the final answer in the required format.
final_answer = solve_rotation_convention()
if final_answer:
    print(f"<<<{final_answer}>>>")
