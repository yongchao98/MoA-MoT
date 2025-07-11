import numpy as np

def Rx(theta_deg):
    """Computes the rotation matrix for a rotation around the X-axis."""
    theta_rad = np.radians(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(theta_deg):
    """Computes the rotation matrix for a rotation around the Y-axis."""
    theta_rad = np.radians(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def Rz(theta_deg):
    """Computes the rotation matrix for a rotation around the Z-axis."""
    theta_rad = np.radians(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# --- 1. Calculate Target Rotation Matrix from Tait-Bryan Angles ---
alpha, beta, gamma = 10, 10, 10
print(f"Step 1: Calculate the target rotation matrix for the extrinsic X-Y-Z convention.")
print(f"Angles: alpha = {alpha} deg, beta = {beta} deg, gamma = {gamma} deg.")
print(f"The calculation is R_target = R_Z(gamma) * R_Y(beta) * R_X(alpha)\n")

R_X_alpha = Rx(alpha)
R_Y_beta = Ry(beta)
R_Z_gamma = Rz(gamma)

R_target = R_Z_gamma @ R_Y_beta @ R_X_alpha

print(f"R_X({alpha}) =\n{R_X_alpha}\n")
print(f"R_Y({beta}) =\n{R_Y_beta}\n")
print(f"R_Z({gamma}) =\n{R_Z_gamma}\n")
print(f"Resulting Target Matrix R_target =\n{R_target}\n")
print("-" * 50)

# --- 2. Test Each Euler Convention ---
alpha_prime, beta_prime, gamma_prime = 139.13, 14.11, -141.05
print(f"Step 2: Find which proper Euler convention with angles:")
print(f"alpha' = {alpha_prime} deg, beta' = {beta_prime} deg, gamma' = {gamma_prime} deg")
print(f"produces the same rotation matrix.\n")

# For an extrinsic convention "ABC" with angles (a,b,c), R = R_C(c) * R_B(b) * R_A(a)
# The angles are associated with the rotation axes in order: a -> A, b -> B, c -> C.
# However, for Euler angles (e.g. ZYZ), the standard is alpha', beta', gamma' corresponding to Z, Y, Z.
# So Z(alpha'), Y(beta'), Z(gamma'). Extrinsic rotation is Rz(gamma') @ Ry(beta') @ Rz(alpha').
conventions_to_test = {
    "A": ("XZX", Rx(gamma_prime) @ Rz(beta_prime) @ Rx(alpha_prime)),
    "C": ("YXY", Ry(gamma_prime) @ Rx(beta_prime) @ Ry(alpha_prime)),
    "D": ("YZY", Ry(gamma_prime) @ Rz(beta_prime) @ Ry(alpha_prime)),
    "E": ("ZYZ", Rz(gamma_prime) @ Ry(beta_prime) @ Rz(alpha_prime)),
    "F": ("ZXZ", Rz(gamma_prime) @ Rx(beta_prime) @ Rz(alpha_prime)),
}

# --- 3. Compare and Find the Match ---
matching_convention_letter = None
for letter, (name, R_test) in conventions_to_test.items():
    # Use a small tolerance for floating point comparison
    if np.allclose(R_target, R_test, atol=1e-4):
        matching_convention_letter = letter
        matching_convention_name = name
        
        print(f"Step 3: A match is found with the '{matching_convention_name}' convention!\n")
        
        # Show the calculation for the matching convention
        if matching_convention_name == 'ZYZ':
             R1 = Rz(alpha_prime)
             R2 = Ry(beta_prime)
             R3 = Rz(gamma_prime)
             print(f"The calculation is R_test = R_Z(gamma') * R_Y(beta') * R_Z(alpha')\n")
             print(f"R_Z({alpha_prime}) =\n{R1}\n")
             print(f"R_Y({beta_prime}) =\n{R2}\n")
             print(f"R_Z({gamma_prime}) =\n{R3}\n")
             print(f"Resulting Test Matrix R_test =\n{R_test}\n")
             
        break

if matching_convention_letter is None:
    print("Step 3: No matching convention was found among the proper Euler choices.")

print("-" * 50)
print(f"The correct convention is '{matching_convention_name}'.")

<<<E>>>