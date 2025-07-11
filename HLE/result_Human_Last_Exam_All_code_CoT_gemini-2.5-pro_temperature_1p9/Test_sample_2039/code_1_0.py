import numpy as np

def Rx(theta):
    """Rotation matrix around X axis (input in radians)"""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(theta):
    """Rotation matrix around Y axis (input in radians)"""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(theta):
    """Rotation matrix around Z axis (input in radians)"""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# --- Input values ---
# Initial extrinsic Tait-Bryan rotation angles in degrees
alpha_deg = 10.0
beta_deg = 10.0
gamma_deg = 10.0

# Equivalent Euler angles in degrees
alpha_p_deg = 139.13
beta_p_deg = 14.11
gamma_p_deg = -141.05

print("Step 1: Formulate the target rotation matrix.")
print("The initial rotation is an extrinsic X-Y-Z rotation (Tait-Bryan angles).")
print("For an extrinsic X-Y-Z rotation, the axes are fixed, and the order of matrix multiplication is reversed.")
print(f"So, the target matrix is R_target = Rz({gamma_deg}°) * Ry({beta_deg}°) * Rx({alpha_deg}°).")

# Convert degrees to radians for numpy functions
alpha_rad = np.deg2rad(alpha_deg)
beta_rad = np.deg2rad(beta_deg)
gamma_rad = np.deg2rad(gamma_deg)

alpha_p_rad = np.deg2rad(alpha_p_deg)
beta_p_rad = np.deg2rad(beta_p_deg)
gamma_p_rad = np.deg2rad(gamma_p_deg)

# Calculate the target rotation matrix
R_target = Rz(gamma_rad) @ Ry(beta_rad) @ Rx(alpha_rad)

print("\nStep 2: Test the candidate Euler angle conventions.")
print("Euler angles represent a sequence of intrinsic rotations, where each rotation is about the new, body-fixed axis.")
print(f"We will test each convention using the given angles: alpha'={alpha_p_deg}°, beta'={beta_p_deg}°, gamma'={gamma_p_deg}°.")

# Define the matrix compositions for each intrinsic convention
# For an intrinsic convention like 'ABC', the matrix is Ra * Rb * Rc
conventions = {
    'A': {'name': 'XZX', 'matrix': Rx(alpha_p_rad) @ Rz(beta_p_rad) @ Rx(gamma_p_rad), 'axes_funcs': ('Rx', 'Rz', 'Rx')},
    'B': {'name': 'XYZ', 'matrix': Rx(alpha_p_rad) @ Ry(beta_p_rad) @ Rz(gamma_p_rad), 'axes_funcs': ('Rx', 'Ry', 'Rz')},
    'C': {'name': 'YXY', 'matrix': Ry(alpha_p_rad) @ Rx(beta_p_rad) @ Ry(gamma_p_rad), 'axes_funcs': ('Ry', 'Rx', 'Ry')},
    'D': {'name': 'YZY', 'matrix': Ry(alpha_p_rad) @ Rz(beta_p_rad) @ Ry(gamma_p_rad), 'axes_funcs': ('Ry', 'Rz', 'Ry')},
    'E': {'name': 'ZYZ', 'matrix': Rz(alpha_p_rad) @ Ry(beta_p_rad) @ Rz(gamma_p_rad), 'axes_funcs': ('Rz', 'Ry', 'Rz')},
    'F': {'name': 'ZXZ', 'matrix': Rz(alpha_p_rad) @ Rx(beta_p_rad) @ Rz(gamma_p_rad), 'axes_funcs': ('Rz', 'Rx', 'Rz')},
}

correct_answer = None
for key, data in conventions.items():
    # Use a tolerance for comparison due to precision of input angles
    if np.allclose(R_target, data['matrix'], atol=1e-4):
        correct_answer = key
        correct_data = data
        break

print("\nStep 3: State the conclusion.")
if correct_answer:
    print(f"The convention that results in an equivalent rotation is '{correct_data['name']}' (Choice {correct_answer}).")
    
    # Construct the final equation string showing the equality
    lhs = f"Rz({gamma_deg}°) * Ry({beta_deg}°) * Rx({alpha_deg}°)"
    
    ax1, ax2, ax3 = correct_data['axes_funcs']
    rhs = f"{ax1}({alpha_p_deg}°) * {ax2}({beta_p_deg}°) * {ax3}({gamma_p_deg}°)"
    
    print("\nThe final equivalence equation is:")
    print(f"{lhs} = {rhs}")
else:
    print("No matching convention was found with the given angles.")

print(f"\n<<<{correct_answer}>>>")