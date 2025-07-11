import numpy as np

def Rx(theta_deg):
    """Creates a 3x3 rotation matrix around the X-axis."""
    theta_rad = np.radians(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(theta_deg):
    """Creates a 3x3 rotation matrix around the Y-axis."""
    theta_rad = np.radians(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(theta_deg):
    """Creates a 3x3 rotation matrix around the Z-axis."""
    theta_rad = np.radians(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

# --- Step 1: Calculate the reference rotation matrix ---
# Extrinsic Tait-Bryan angles for R_ref: X_alpha Y_beta Z_gamma
# The composite matrix is R_z(gamma) @ R_y(beta) @ R_x(alpha)
alpha, beta, gamma = 10.0, 10.0, 10.0

R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# --- Step 2: Define candidate conventions and angles ---
# Given Euler angles for the candidate conventions
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# --- Step 3: Test each convention ---
# These are intrinsic rotations, so matrix multiplication is in order of convention
# e.g., for ZYZ: R_z(alpha') @ R_y(beta') @ R_z(gamma')
candidates = {
    'A. XZX': Rx(alpha_p) @ Rz(beta_p) @ Rx(gamma_p),
    'B. XYZ': Rx(alpha_p) @ Ry(beta_p) @ Rz(gamma_p),
    'C. YXY': Ry(alpha_p) @ Rx(beta_p) @ Ry(gamma_p),
    'D. YZY': Ry(alpha_p) @ Rz(beta_p) @ Ry(gamma_p),
    'E. ZYZ': Rz(alpha_p) @ Ry(beta_p) @ Rz(gamma_p),
    'F. ZXZ': Rz(alpha_p) @ Rx(beta_p) @ Rz(gamma_p)
}

# --- Step 4: Compare and find the match ---
matching_convention = None
for name, R_cand in candidates.items():
    # Use np.allclose to compare matrices with a tolerance for floating point errors
    if np.allclose(R_ref, R_cand, atol=1e-3):
        matching_convention = name
        break

# --- Step 5: Output the result ---
if matching_convention:
    print(f"A match was found for convention: {matching_convention}")
    print("\nThe original rotation is given by the extrinsic Tait-Bryan angles:")
    print(f"α = {alpha}°, β = {beta}°, γ = {gamma}° following the X-Y-Z convention.")
    print("This results in the composite rotation matrix R_ref = Rz(γ) @ Ry(β) @ Rx(α).\n")

    print("The equivalent rotation is described by the proper Euler angles:")
    print(f"α' = {alpha_p}°, β' = {beta_p}°, γ' = {gamma_p}° using the {matching_convention.split('.')[1].strip()} convention.")
    print(f"This results in the composite rotation matrix R_{matching_convention.split('.')[1].strip()} = Rz(α') @ Ry(β') @ Rz(γ').\n")

    print("The final check shows that the two matrices are equivalent:")
    print("Reference Matrix (from X_10 Y_10 Z_10):")
    print(R_ref)
    print(f"\nMatching Candidate Matrix ({matching_convention}):")
    print(candidates[matching_convention])
else:
    print("No matching convention was found among the options.")
