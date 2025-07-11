import cvxpy as cp
import numpy as np
import itertools

def solve_for_z(A):
    """
    Calculates the smallest z for a given matrix A such that zB - A is PSD
    for some "nice" matrix B.
    """
    n = A.shape[0]
    
    # 1. Generate all 2^n configurations of {-1, 1}^n vectors
    configs = list(itertools.product([-1, 1], repeat=n))
    num_configs = len(configs)
    
    # 2. Define optimization variables
    # z is the scalar we want to minimize
    z = cp.Variable()
    # p is the probability distribution over the 2^n configurations
    p = cp.Variable(num_configs, nonneg=True)

    # 3. Construct the "nice" matrix B
    # B is a convex combination of outer products xx^T
    B = sum(p[i] * np.outer(configs[i], configs[i]) for i in range(num_configs))
    
    # 4. Define the constraints
    constraints = [
        # p must be a probability distribution
        cp.sum(p) == 1,
        # The random variables must be unbiased, i.e., E[X_i] = 0 for all i.
        # This is equivalent to E[X] = 0, where X is the random vector.
        sum(p[i] * np.array(configs[i]) for i in range(num_configs)) == np.zeros(n),
        # The core constraint: zB - A must be positive semidefinite
        z * B - A >> 0
    ]
    
    # 5. Define and solve the problem
    objective = cp.Minimize(z)
    problem = cp.Problem(objective, constraints)
    problem.solve()

    # 6. Print the results
    if problem.status == cp.OPTIMAL:
        z_val = z.value
        B_val = B.value
        C_val = z_val * B_val - A
        
        np.set_printoptions(precision=4, suppress=True)
        print("For the given matrix A:")
        print(A)
        print("\nThe smallest value of z found is: {:.4f}\n".format(z_val))
        
        print("The final equation is A = z*B - C\n")

        print("z = \n{:.4f}\n".format(z_val))
        
        print("B (nice matrix) =")
        print(B_val)
        print("\nC (positive semidefinite matrix) =")
        print(C_val)
        
        # Verify C is PSD by checking its eigenvalues
        eigenvalues = np.linalg.eigvalsh(C_val)
        print("\nEigenvalues of C (should be non-negative):")
        print(eigenvalues)

    else:
        print("The problem could not be solved to optimality.")
        print("Status:", problem.status)

# Define a 3x3 correlation matrix A that is not "nice"
n = 3
A = np.full((n, n), 0.5)
np.fill_diagonal(A, 1)

solve_for_z(A)