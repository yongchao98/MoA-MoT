import numpy as np

def Rx(degrees):
    """Computes the 3x3 rotation matrix for a rotation around the X-axis."""
    rad = np.deg2rad(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[1, 0,  0],
                     [0, c, -s],
                     [0, s,  c]])

def Ry(degrees):
    """Computes the 3x3 rotation matrix for a rotation around the Y-axis."""
    rad = np.deg2rad(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(degrees):
    """Computes the 3x3 rotation matrix for a rotation around the Z-axis."""
    rad = np.deg2rad(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, -s, 0],
                     [s,  c, 0],
                     [0,  0, 1]])

# 1. Define the initial Tait-Bryan rotation angles.
alpha_tb = 10.0
beta_tb = 10.0
gamma_tb = 10.0

# 2. Calculate the reference rotation matrix using the X-Y-Z convention.
# For an extrinsic XYZ rotation with angles a, b, g, the matrix is R_x(a)R_y(b)R_z(g).
R_ref = Rx(alpha_tb) @ Ry(beta_tb) @ Rz(gamma_tb)

print("Step 1: Calculate the rotation matrix from the initial Tait-Bryan XYZ rotation.")
print(f"R_ref = Rx({alpha_tb}) @ Ry({beta_tb}) @ Rz({gamma_tb})")
print("Reference Matrix (from Tait-Bryan angles):")
print(R_ref)
print("-" * 50)

# 3. Define the candidate proper Euler angles.
alpha_e = 139.13
beta_e = 14.11
gamma_e = -141.05

print("Step 2: Test candidate Euler angle conventions with the following angles:")
print(f"alpha' = {alpha_e}, beta' = {beta_e}, gamma' = {gamma_e}")
print("-" * 50)


# 4. Define the conventions to test (name, axis sequence).
# For an extrinsic rotation with convention ABC and angles a,b,c the matrix is Ra(a)Rb(b)Rc(c)
conventions_to_test = {
    "A. XZX": (Rx, Rz, Rx),
    "C. YXY": (Ry, Rx, Ry),
    "D. YZY": (Ry, Rz, Ry),
    "E. ZYZ": (Rz, Ry, Rz),
    "F. ZXZ": (Rz, Rx, Rz),
}

# 5. Iterate, calculate, and compare.
found_match = False
for name, (R1, R2, R3) in conventions_to_test.items():
    
    # Construct the final equation string for printing
    r_map = {Rx: 'Rx', Ry: 'Ry', Rz: 'Rz'}
    eq_str = f"R_cand = {r_map[R1]}({alpha_e}) @ {r_map[R2]}({beta_e}) @ {r_map[R3]}({gamma_e})"

    # Calculate the candidate matrix
    R_candidate = R1(alpha_e) @ R2(beta_e) @ R3(gamma_e)
    
    # Check if the candidate matrix is close to the reference matrix
    # Using a small tolerance to account for floating point and rounding of input angles.
    is_match = np.allclose(R_ref, R_candidate, atol=1e-3)
    
    print(f"Testing convention: {name}")
    print(eq_str)
    print("Resulting Matrix:")
    print(R_candidate)

    if is_match:
        print("\n>>> SUCCESS: This matrix matches the reference matrix.")
        print(f"The correct convention is {name}.")
        found_match = True
        break
    else:
        print(">>> This matrix does not match the reference matrix.\n")
        
if not found_match:
    print("No matching convention was found among the provided proper Euler choices.")
