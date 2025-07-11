import numpy as np
import cvxpy as cp

def solve_for_z(A):
    """
    Finds the minimum z for a given matrix A such that A = z*B - C,
    where B is in the 3x3 elliptope and C is positive semidefinite.
    """
    print(f"Finding the decomposition for matrix A:\n{A}\n")
    
    # B is the variable we are solving for, a 3x3 symmetric matrix.
    B_var = cp.Variable((3, 3), symmetric=True)
    # z_param is a parameter we will set during a binary search.
    z_param = cp.Parameter(nonneg=True)

    # Constraints for B to be in the 3x3 elliptope (E_3).
    # This means B must be a correlation matrix satisfying the triangle inequalities.
    constraints = [
        B_var >> 0,  # Must be positive semidefinite
        cp.diag(B_var) == 1,  # Must have 1s on the diagonal
        1 + B_var[0,1] + B_var[0,2] + B_var[1,2] >= 0,
        1 + B_var[0,1] - B_var[0,2] - B_var[1,2] >= 0,
        1 - B_var[0,1] + B_var[0,2] - B_var[1,2] >= 0,
        1 - B_var[0,1] - B_var[0,2] + B_var[1,2] >= 0,
    ]
    # The main constraint: z*B - A must be positive semidefinite.
    constraints.append(z_param * B_var - A >> 0)

    # Create the feasibility problem. The objective is irrelevant.
    prob = cp.Problem(cp.Minimize(0), constraints)

    # Binary search for the minimal feasible z.
    # K_G is known to be between 1 and 2.
    low = 1.0
    high = 2.0 
    print("Searching for the minimum z for this specific matrix...")
    # Iterate to find z with high precision.
    for _ in range(100):
        mid = (low + high) / 2
        z_param.value = mid
        # solve() returns the objective value, or inf if infeasible.
        prob.solve(solver=cp.SCS, warm_start=True)
        if prob.status in ['optimal', 'optimal_inaccurate']:
            high = mid
        else:
            low = mid
            
    z_min = high
    print(f"Found minimum z: {z_min:.8f}\n")

    # Resolve with the found z_min to get the matrices B and C.
    z_param.value = z_min
    prob.solve(solver=cp.SCS)
    B_sol = B_var.value
    C_sol = z_min * B_sol - A
    
    # The final equation is A = z * B - C. Let's print all numbers.
    print("The final equation A = z * B - C with computed values:")
    print("-" * 55)
    # Using np.array2string for nice formatting
    print(f"A = \n{np.array2string(A, precision=4)}\n")
    print(f"z = {z_min:.8f}\n")
    print(f"B = \n{np.array2string(B_sol, precision=4)}\n")
    print(f"C = z*B - A = \n{np.array2string(C_sol, precision=4)}\n")
    print("-" * 55)

    # Verify that the resulting C matrix is indeed positive semidefinite.
    eigenvalues_C = np.linalg.eigvalsh(C_sol)
    print(f"Verification: Eigenvalues of C are {np.array2string(eigenvalues_C, precision=4)}")
    # Small negative eigenvalues can occur due to numerical tolerance.
    if np.all(eigenvalues_C >= -1e-7):
        print("This confirms C is positive semidefinite.")
    else:
        print("Warning: C is not positive semidefinite, there might be a numerical issue.")

if __name__ == '__main__':
    # We choose a specific 3x3 correlation matrix A that is not "nice"
    # (i.e., not in the elliptope). The computed z will be a lower bound for K_G.
    # For A to not be in E_3, it must violate one of the triangle inequalities.
    # For A[0,1]=0.4, A[0,2]=0.4, A[1,2]=-0.6: 1-0.4-0.4+(-0.6) = -0.4 < 0.
    A_matrix = np.array([
        [1.0, 0.4, 0.4],
        [0.4, 1.0, -0.6],
        [0.4, -0.6, 1.0]
    ])
    solve_for_z(A_matrix)