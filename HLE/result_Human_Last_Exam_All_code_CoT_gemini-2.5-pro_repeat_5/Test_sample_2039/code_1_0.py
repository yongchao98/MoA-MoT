import numpy as np
from scipy.spatial.transform import Rotation

def get_rot_matrix_string(rot_matrix):
    """Formats a 3x3 numpy array into a readable string."""
    return '\n'.join(['[' + ' '.join([f'{cell:8.4f}' for cell in row]) + ']' for row in rot_matrix])

# --- Step 1: Calculate the Target Rotation Matrix ---
# Extrinsic rotation X_alpha Y_beta Z_gamma with alpha=beta=gamma=10 deg.
# This means R = R_Z(gamma) * R_Y(beta) * R_X(alpha).
# In SciPy's 'XYZ' extrinsic notation, this corresponds to angles=[alpha, beta, gamma].
alpha, beta, gamma = 10.0, 10.0, 10.0
target_angles_deg = [alpha, beta, gamma]
target_convention_str = "extrinsic XYZ"

# Create the rotation object
r_target = Rotation.from_euler('XYZ', target_angles_deg, degrees=True)
R_target = r_target.as_matrix()

print("Step 1: Calculate the target rotation matrix from the given Tait-Bryan angles.")
print(f"The convention is {target_convention_str} with angles:")
print(f"alpha = {alpha:.2f} deg (around X)")
print(f"beta  = {beta:.2f} deg (around Y)")
print(f"gamma = {gamma:.2f} deg (around Z)")
print(f"\nThe corresponding rotation equation is R_target = Rz({gamma:.2f}) * Ry({beta:.2f}) * Rx({alpha:.2f})")
print("Resulting R_target matrix:")
print(get_rot_matrix_string(R_target))
print("-" * 40)


# --- Step 2: Define the Euler Angles and Conventions to Test ---
euler_angles_deg = [139.13, 14.11, -141.05]
alpha_p, beta_p, gamma_p = euler_angles_deg

# The answer choices are the conventions to test.
# "Proper Euler angles" implies an intrinsic rotation sequence (e.g., z-y'-z'').
conventions_to_test = {
    "A": "XZX", "B": "XYZ", "C": "YXY",
    "D": "YZY", "E": "ZYZ", "F": "ZXZ"
}

print("Step 2: Test the given proper Euler angles with different intrinsic conventions.")
print("The angles to test are:")
print(f"alpha' = {alpha_p:.2f} deg")
print(f"beta'  = {beta_p:.2f} deg")
print(f"gamma' = {gamma_p:.2f} deg")
print("-" * 40)


# --- Step 3 & 4: Iterate, Calculate, and Compare ---
print("Step 3: Finding the matching convention...\n")
final_answer_key = None

for key, conv_name in conventions_to_test.items():
    # Assume intrinsic rotation (lowercase in scipy)
    r_test = Rotation.from_euler(conv_name.lower(), euler_angles_deg, degrees=True)
    R_test = r_test.as_matrix()

    # Check if the matrices are close enough
    is_match = np.allclose(R_target, R_test, atol=1e-4)

    print(f"--- Testing Convention {key}: {conv_name} (intrinsic) ---")
    axis1, axis2, axis3 = conv_name[0], conv_name[1], conv_name[2]
    # For intrinsic rotation R = R_axis1(alpha') * R_axis2(beta') * R_axis3(gamma')
    print(f"The rotation equation is R_test = R{axis1.lower()}({alpha_p:.2f}) * R{axis2.lower()}({beta_p:.2f}) * R{axis3.lower()}({gamma_p:.2f})")
    print("Resulting R_test matrix:")
    print(get_rot_matrix_string(R_test))

    if is_match:
        print("\n>>> MATCH FOUND! <<<\n")
        print(f"The intrinsic {conv_name} convention produces the same rotation matrix.")
        final_answer_key = key
        break 
    else:
        print("...No match.\n")

if final_answer_key is None:
    print("No matching convention was found among the options.")

<<<E>>>