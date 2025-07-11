import numpy as np
import cvxpy as cp

def solve_for_A(A):
    """
    Solves the SDP to find the smallest z for a given correlation matrix A.
    min z  s.t. A = zB - C, B in CUT_n, C >= 0
    This is equivalent to the SDP:
    min sum(q) s.t. sum(q_k * B_k) - A >= 0, q_k >= 0
    where z = sum(q) and B = (1/z) * sum(q_k * B_k)
    """
    n = A.shape[0]
    
    # Generate the basis matrices for CUT_n
    # There are 2^(n-1) such matrices
    cut_basis = []
    for i in range(2**(n - 1)):
        x = np.ones(n)
        # Create a vector x from {-1, 1}^n
        for j in range(n - 1):
            if (i >> j) & 1:
                x[j + 1] = -1
        B_k = np.outer(x, x)
        cut_basis.append(B_k)

    # Define the SDP variables
    q = cp.Variable(len(cut_basis), nonneg=True)
    
    # Form the matrix B_tilde = sum(q_k * B_k)
    B_tilde = sum(q[k] * cut_basis[k] for k in range(len(cut_basis)))
    
    # Define the constraints
    constraints = [B_tilde - A >> 0]
    
    # Define the objective function
    objective = cp.Minimize(cp.sum(q))
    
    # Define and solve the problem
    problem = cp.Problem(objective, constraints)
    problem.solve()
    
    # Extract results
    z_val = problem.value
    q_val = q.value
    
    # Reconstruct B and C
    B_val = sum(q_val[k] * cut_basis[k] for k in range(len(cut_basis))) / z_val
    C_val = z_val * B_val - A
    
    print(f"For the given matrix A:\n{A}\n")
    print(f"The smallest z is: {z_val:.4f}\n")
    
    print("This z is a lower bound for K_G.")
    print("The equation is A = z*B - C\n")
    
    print(f"z = {z_val:.4f}\n")
    
    print(f"B (nice matrix):\n{B_val}\n")
    
    print(f"C (positive semidefinite matrix):\n{C_val}\n")
    
    # Verification
    print("Verification: z*B - C\n", z_val * B_val - C_val)
    print("\nEigenvalues of C (should be non-negative):")
    eigvals_C = np.linalg.eigvalsh(C_val)
    # Set very small values to 0 for cleaner output
    eigvals_C[np.isclose(eigvals_C, 0)] = 0
    print(eigvals_C)


# Create a 3x3 correlation matrix A from vectors at 120 degrees
A = np.array([
    [1.0, -0.5, -0.5],
    [-0.5, 1.0, -0.5],
    [-0.5, -0.5, 1.0]
])

solve_for_A(A)
