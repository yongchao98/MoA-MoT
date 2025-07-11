import numpy as np
from scipy.spatial.transform import Rotation

# Step 1: Calculate the rotation matrix for the Tait-Bryan angles
# Extrinsic X, Y, Z with angles 10, 10, 10 degrees.
# Scipy's 'xyz' is intrinsic, so to get extrinsic XYZ, we use intrinsic 'zyx' in reverse order of angles.
# However, it's more direct to construct it from basic rotation matrices to avoid confusion.
# Extrinsic R = Rz(10) * Ry(10) * Rx(10)
r_x = Rotation.from_euler('x', 10, degrees=True).as_matrix()
r_y = Rotation.from_euler('y', 10, degrees=True).as_matrix()
r_z = Rotation.from_euler('z', 10, degrees=True).as_matrix()
r_tb_matrix = r_z @ r_y @ r_x

# Step 2: Define the target Euler angles
alpha_prime = 139.13
beta_prime = 14.11
gamma_prime = -141.05

# Step 3: Test each proper Euler angle convention
conventions = {
    "A": "xzx",
    "B": "xyz", # This is Tait-Bryan, not proper
    "C": "yxy",
    "D": "yzy",
    "E": "zyz",
    "F": "zxz",
}

results = {}

for key, conv in conventions.items():
    if key == "B": continue # Skip the non-proper convention
    # Scipy conventions are intrinsic.
    r = Rotation.from_euler(conv, [alpha_prime, beta_prime, gamma_prime], degrees=True)
    results[key] = r.as_matrix()

# Step 4: Compare the matrices
print("Reference Tait-Bryan Matrix (X-Y-Z, 10-10-10):")
print(r_tb_matrix)
print("\n" + "="*30 + "\n")

for key, matrix in results.items():
    print(f"Matrix for convention {conventions[key].upper()} with angles {alpha_prime, beta_prime, gamma_prime}:")
    print(matrix)
    # Check if the matrices are close
    is_close = np.allclose(r_tb_matrix, matrix)
    print(f"Is it equivalent? {is_close}")
    print("\n" + "="*30 + "\n")
