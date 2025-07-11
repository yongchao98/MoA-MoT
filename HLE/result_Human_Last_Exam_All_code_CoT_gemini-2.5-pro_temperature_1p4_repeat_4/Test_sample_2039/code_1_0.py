import numpy as np

def Rx(deg):
    """Computes the rotation matrix around the X-axis for a given angle in degrees."""
    rad = np.deg2rad(deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(deg):
    """Computes the rotation matrix around the Y-axis for a given angle in degrees."""
    rad = np.deg2rad(deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(deg):
    """Computes the rotation matrix around the Z-axis for a given angle in degrees."""
    rad = np.deg2rad(deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# --- Step 1: Define the initial Tait-Bryan rotation and calculate its matrix ---
alpha_tb, beta_tb, gamma_tb = 10.0, 10.0, 10.0
tb_convention = "XYZ"

print(f"Calculating the reference rotation matrix for the extrinsic Tait-Bryan convention {tb_convention}...")
print(f"The rotation angles are: alpha = {alpha_tb} deg, beta = {beta_tb} deg, gamma = {gamma_tb} deg.")
# For extrinsic XYZ, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(gamma_tb) @ Ry(beta_tb) @ Rx(alpha_tb)
print("\nReference Rotation Matrix (Tait-Bryan XYZ):")
print(R_ref)

# --- Step 2: Define the target Euler angles and the possible conventions ---
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05
print(f"\nNow, testing proper Euler angle conventions with the angles:")
print(f"alpha' = {alpha_p} deg, beta' = {beta_p} deg, gamma' = {gamma_p} deg.")

conventions = {
    'A': 'XZX',
    'B': 'XYZ', # Included for completeness, though it's Tait-Bryan
    'C': 'YXY',
    'D': 'YZY',
    'E': 'ZYZ',
    'F': 'ZXZ'
}
rot_funcs = {'X': Rx, 'Y': Ry, 'Z': Rz}
found_match = False
correct_answer = None

# --- Steps 3 & 4: Iterate through conventions, calculate matrices, and compare ---
for key, conv_str in conventions.items():
    if len(set(conv_str)) == 3 and conv_str[0] == conv_str[2]:
        print(f"\n--- This case is impossible for proper Euler Angles. {conv_str}")

    # Get the rotation functions for the current convention
    R_func1 = rot_funcs[conv_str[0]]
    R_func2 = rot_funcs[conv_str[1]]
    R_func3 = rot_funcs[conv_str[2]]

    # Calculate the candidate matrix for the extrinsic convention
    # R = R_axis3(gamma') * R_axis2(beta') * R_axis1(alpha')
    R_cand = R_func3(gamma_p) @ R_func2(beta_p) @ R_func1(alpha_p)

    # Compare the candidate matrix with the reference matrix
    if np.allclose(R_ref, R_cand):
        found_match = True
        correct_answer = key
        print(f"\nSUCCESS: A match is found for convention '{conv_str}' (Choice {key}).")
        print("The two rotation sequences produce the same orientation.")
        print("\nThe final equivalent rotation equation is:")
        
        eq_str1 = f"R_{tb_convention[2]}({gamma_tb}) * R_{tb_convention[1]}({beta_tb}) * R_{tb_convention[0]}({alpha_tb})"
        eq_str2 = f"R_{conv_str[2]}({gamma_p}) * R_{conv_str[1]}({beta_p}) * R_{conv_str[0]}({alpha_p})"
        
        print(f"{eq_str1}  ==  {eq_str2}")
        break

if not found_match:
    print("\nCould not find a matching convention among the options.")

print(f"\nBased on the calculation, the correct convention is {conventions.get(correct_answer)}.")
print("<<<F>>>")