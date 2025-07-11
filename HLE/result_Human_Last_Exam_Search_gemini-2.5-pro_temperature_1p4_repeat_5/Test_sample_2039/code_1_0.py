import numpy as np
from scipy.spatial.transform import Rotation as R

# Define the initial Tait-Bryan rotation
# Extrinsic X-Y-Z rotation with alpha=beta=gamma=10 degrees
alpha_tb = 10
beta_tb = 10
gamma_tb = 10
# For scipy, extrinsic 'xyz' means the order of operations is Rz, then Ry, then Rx
# The problem states X_alpha Y_beta Z_gamma which corresponds to R = Rx(a) * Ry(b) * Rz(g)
# This corresponds to an intrinsic 'zyx' rotation in scipy's convention or an extrinsic 'xyz' applied in order z, y, x
# Let's build it manually to be sure.
# Rx(10) * Ry(10) * Rz(10)
r_initial = R.from_euler('x', 10, degrees=True) * R.from_euler('y', 10, degrees=True) * R.from_euler('z', 10, degrees=True)
initial_matrix = r_initial.as_matrix()

# Define the equivalent proper Euler angles
alpha_p = 139.13
beta_p = 14.11
gamma_p = -141.05

# Define the candidate conventions
conventions = {
    "A": "xzx",
    "C": "yxy",
    "D": "yzy",
    "E": "zyz",
    "F": "zxz"
}

# Calculate the rotation matrix for each convention and compare
results = {}
for key, seq in conventions.items():
    # Scipy's 'euler' method uses intrinsic rotations. The order of letters is the order of rotations.
    # e.g., 'zyz' means rotate by alpha' around Z, then by beta' around the new Y', then by gamma' around the new Z''.
    r_candidate = R.from_euler(seq, [alpha_p, beta_p, gamma_p], degrees=True)
    candidate_matrix = r_candidate.as_matrix()
    # Compare the candidate matrix with the initial matrix
    is_close = np.allclose(initial_matrix, candidate_matrix)
    results[key] = is_close

print("Initial Tait-Bryan (X_10 Y_10 Z_10) Matrix:")
print(initial_matrix)
print("\nComparison Results:")
print(results)
