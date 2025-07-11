import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_problem():
    """
    Finds the Euler angle convention that matches a given extrinsic Tait-Bryan rotation.
    """
    # 1. Define the initial extrinsic Tait-Bryan rotation
    initial_angles_deg = [10, 10, 10]
    # The problem specifies extrinsic rotation using X, then Y, then Z axis.
    # In scipy, for extrinsic rotation, this corresponds to the sequence 'xyz'.
    # The resulting matrix is Rz(gamma) * Ry(beta) * Rx(alpha).
    initial_rotation = Rotation.from_euler('xyz', initial_angles_deg, degrees=True, extrinsic=True)
    initial_matrix = initial_rotation.as_matrix()

    print(f"The initial rotation is an extrinsic XYZ rotation with angles:\n"
          f"alpha = {initial_angles_deg[0]}°, beta = {initial_angles_deg[1]}°, gamma = {initial_angles_deg[2]}°\n")

    # 2. Define the target Euler angles and the candidate conventions
    target_angles_deg = [139.13, 14.11, -141.05]
    conventions = {
        'A': 'xzx',
        'B': 'xyz',
        'C': 'yxy',
        'D': 'yzy',
        'E': 'zyz',
        'F': 'zxz'
    }

    print(f"We are searching for an equivalent intrinsic Euler convention with angles:\n"
          f"alpha' = {target_angles_deg[0]}°, beta' = {target_angles_deg[1]}°, gamma' = {target_angles_deg[2]}°\n")

    # 3. Iterate through conventions, create rotation matrices, and compare them
    found_convention_key = None
    found_convention_str = ""

    print("Testing candidate conventions:")
    for key, conv_str in conventions.items():
        # Standard Euler angles are intrinsic rotations. This is the default in scipy.
        candidate_rotation = Rotation.from_euler(conv_str, target_angles_deg, degrees=True)
        candidate_matrix = candidate_rotation.as_matrix()

        # Compare the matrices with a tolerance to account for precision of input angles
        if np.allclose(initial_matrix, candidate_matrix, atol=1e-4):
            found_convention_key = key
            found_convention_str = conv_str
            print(f"  - Convention {key} ({conv_str.upper()}): Match found!")
            break
        else:
            print(f"  - Convention {key} ({conv_str.upper()}): No match.")

    # 4. Output the result
    if found_convention_key:
        print("\n---")
        print("Conclusion: The correct convention is ZYZ.")
        print("The final equivalence is:")
        print(
            f"Extrinsic_XYZ({initial_angles_deg[0]:.2f}°, {initial_angles_deg[1]:.2f}°, {initial_angles_deg[2]:.2f}°) == "
            f"Intrinsic_{found_convention_str.upper()}({target_angles_deg[0]:.2f}°, {target_angles_deg[1]:.2f}°, {target_angles_deg[2]:.2f}°)"
        )
        print("---")
        # Returning the answer key
        return found_convention_key
    else:
        print("\nNo matching convention was found among the given options.")
        return None

# Execute the solver
result_key = solve_rotation_problem()
print(f"\n<<<E>>>")