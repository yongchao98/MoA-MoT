import numpy as np

def Rx(deg):
    """Rotation matrix around X-axis for an angle in degrees."""
    rad = np.radians(deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(deg):
    """Rotation matrix around Y-axis for an angle in degrees."""
    rad = np.radians(deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def Rz(deg):
    """Rotation matrix around Z-axis for an angle in degrees."""
    rad = np.radians(deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# --- Step 1: Calculate the reference rotation matrix ---
alpha_deg = 10
beta_deg = 10
gamma_deg = 10

# For extrinsic X-Y-Z, the matrix is R = Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(gamma_deg) @ Ry(beta_deg) @ Rx(alpha_deg)

print("The initial rotation is given by extrinsic Tait-Bryan angles:")
print(f"alpha = {alpha_deg}°, beta = {beta_deg}°, gamma = {gamma_deg}° in X-Y-Z order.")
print("\nThe corresponding reference rotation matrix R_ref = R_z(10) @ R_y(10) @ R_x(10) is:")
print(R_ref)
print("\n" + "-"*50)


# --- Step 2: Test each Euler angle convention ---
alpha_p_deg = 139.13
beta_p_deg = 14.11
gamma_p_deg = -141.05

print("Now, we test which convention with the angles:")
print(f"alpha' = {alpha_p_deg}°, beta' = {beta_p_deg}°, gamma' = {gamma_p_deg}°")
print("produces the same rotation matrix.")

# Define the conventions and their matrix calculations
conventions = {
    "A": ("XZX", Rx(alpha_p_deg) @ Rz(beta_p_deg) @ Rx(gamma_p_deg)),
    "B": ("XYZ", Rx(alpha_p_deg) @ Ry(beta_p_deg) @ Rz(gamma_p_deg)),
    "C": ("YXY", Ry(alpha_p_deg) @ Rx(beta_p_deg) @ Ry(gamma_p_deg)),
    "D": ("YZY", Ry(alpha_p_deg) @ Rz(beta_p_deg) @ Ry(gamma_p_deg)),
    "E": ("ZYZ", Rz(alpha_p_deg) @ Ry(beta_p_deg) @ Rz(gamma_p_deg)),
    "F": ("ZXZ", Rz(alpha_p_deg) @ Rx(beta_p_deg) @ Rz(gamma_p_deg)),
}

# --- Step 3: Compare and find the match ---
correct_option = None
for option, (name, R_test) in conventions.items():
    # Check if the calculated matrix is close to the reference matrix
    if np.allclose(R_ref, R_test, atol=1e-4):
        correct_option = option
        print(f"\nMatch found! The correct convention is {option}. {name}.")
        print(f"\nThe equation for the {name} convention is:")
        print(f"R_test = R_{name[0]}({alpha_p_deg}°) @ R_{name[1]}({beta_p_deg}°) @ R_{name[2]}({gamma_p_deg}°)")
        print("\nThis produces the matrix:")
        print(R_test)
        print("\nThis matrix is numerically equivalent to R_ref, confirming it is the correct convention.")
        break

if correct_option is None:
    print("\nCould not find a matching convention among the options.")
