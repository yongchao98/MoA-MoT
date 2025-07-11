import numpy as np

def Rx(theta_deg):
    """Computes the rotation matrix for a rotation around the X-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(theta_deg):
    """Computes the rotation matrix for a rotation around the Y-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def Rz(theta_deg):
    """Computes the rotation matrix for a rotation around the Z-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# --- Step 1: Calculate the reference rotation matrix ---
alpha, beta, gamma = 10, 10, 10
print(f"Calculating the reference rotation matrix for extrinsic X({alpha}) Y({beta}) Z({gamma}) rotation.")
# For extrinsic XYZ, the matrix is Rz(gamma) @ Ry(beta) @ Rx(alpha)
R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)
print("Reference Matrix:\n", R_ref)
print("-" * 40)

# --- Step 2 & 3: Test each convention and find the match ---
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05
angles_p = (alpha_p, beta_p, gamma_p)
print(f"Testing Euler angle conventions with angles α'={alpha_p}°, β'={beta_p}°, γ'={gamma_p}°\n")

# Define the conventions to test (assuming intrinsic rotations)
conventions = {
    "XZX": lambda a, b, c: Rx(a) @ Rz(b) @ Rx(c),
    "YXY": lambda a, b, c: Ry(a) @ Rx(b) @ Ry(c),
    "YZY": lambda a, b, c: Ry(a) @ Rz(b) @ Ry(c),
    "ZYZ": lambda a, b, c: Rz(a) @ Ry(b) @ Rz(c),
    "ZXZ": lambda a, b, c: Rz(a) @ Rx(b) @ Rz(c),
}

# We also check the Tait-Bryan XYZ convention for completeness
# For intrinsic XYZ, the matrix is Rx(a) @ Ry(b) @ Rz(c)
all_choices = conventions.copy()
all_choices["XYZ"] = lambda a, b, c: Rx(a) @ Ry(b) @ Rz(c)

found_match = False
for name, func in all_choices.items():
    print(f"Testing convention: {name}")
    R_test = func(*angles_p)
    print(f"Resulting Matrix for {name}:\n", R_test)
    
    # Compare with reference matrix using a tolerance for floating point errors
    # The tolerance is set based on the precision of the input angles
    if np.allclose(R_ref, R_test, atol=1e-4):
        print(f"\nSUCCESS: The matrix for the {name} convention matches the reference matrix.")
        print(f"The correct Euler angle convention is {name}.")
        found_match = True
        break
    else:
        print("Matrices do not match.\n")

if not found_match:
    print("\nCould not find a matching convention among the choices.")
