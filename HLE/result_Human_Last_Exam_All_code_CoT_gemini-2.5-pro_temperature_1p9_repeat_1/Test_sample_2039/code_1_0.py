import numpy as np

def Rx(theta_deg):
    """Rotation matrix around X-axis"""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(theta_deg):
    """Rotation matrix around Y-axis"""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(theta_deg):
    """Rotation matrix around Z-axis"""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

# --- Step 1: Calculate the reference rotation matrix from Tait-Bryan angles ---
alpha = 10  # degrees
beta = 10   # degrees
gamma = 10  # degrees

# For extrinsic XYZ rotation, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# --- Step 2: Define the Euler angles for testing ---
alpha_p = 139.13  # degrees
beta_p = 14.11    # degrees
gamma_p = -141.05 # degrees

# --- Step 3 & 4: Test each proper Euler convention and compare ---
conventions = {
    "A": ("XZX", Rx(gamma_p) @ Rz(beta_p) @ Rx(alpha_p)),
    "C": ("YXY", Ry(gamma_p) @ Rx(beta_p) @ Ry(alpha_p)),
    "D": ("YZY", Ry(gamma_p) @ Rz(beta_p) @ Ry(alpha_p)),
    "E": ("ZYZ", Rz(gamma_p) @ Ry(beta_p) @ Rz(alpha_p)),
    "F": ("ZXZ", Rz(gamma_p) @ Rx(beta_p) @ Rz(alpha_p)),
}

found = False
for choice, (name, R_test) in conventions.items():
    # np.allclose checks if two arrays are element-wise equal within a tolerance
    if np.allclose(R_ref, R_test):
        print(f"Match found for convention: {name}")
        print("This corresponds to an extrinsic rotation R = R_{axis3}(gamma') R_{axis2}(beta') R_{axis1}(alpha')")
        print("The final equation is formed by the convention:")
        print(f"{name[2]}({gamma_p}°) {name[1]}({beta_p}°) {name[0]}({alpha_p}°)")
        print(f"\nThe corresponding answer choice is {choice}.")
        found = True
        break

if not found:
    print("No matching convention was found among the given options.")

<<<E>>>