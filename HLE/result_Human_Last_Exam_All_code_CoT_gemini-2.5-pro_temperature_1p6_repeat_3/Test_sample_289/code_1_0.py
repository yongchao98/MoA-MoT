import numpy as np

def solve_eigenvalue_problem():
    """
    This function demonstrates the solution to the eigenvalue problem.

    The problem asks for the largest size of a set S of non-real eigenvalues
    for a matrix A where A^3 = A*.
    The derivation shows that any such eigenvalue lambda must satisfy lambda^3 = conj(lambda).
    The non-real solutions to this equation are i and -i.

    This script verifies that a matrix with these eigenvalues satisfies the condition
    and therefore, the largest possible size of S is 2.
    """
    print("Step 1: The non-real eigenvalues must be i and -i, as derived from lambda^3 = conj(lambda).")
    
    # The set of potential non-real eigenvalues
    non_real_solutions = [1j, -1j]
    
    print("\nStep 2: Construct a matrix A with these non-real eigenvalues.")
    # For simplicity, we can construct a diagonal matrix. Any n >= 2 works.
    # We can also add any of the real solutions (0, 1, -1) to the diagonal.
    # Let's use n=4 and include all possible eigenvalues.
    all_possible_eigenvalues = [1j, -1j, 1, -1, 0]
    # We can only have n eigenvalues for an n x n matrix.
    eigenvalues = np.array([1j, -1j, 1, 0])
    A = np.diag(eigenvalues)
    
    print("Constructed a 4x4 matrix A with eigenvalues:", eigenvalues)
    print("A =")
    print(A)

    print("\nStep 3: Verify that the constructed matrix A satisfies A^3 = A*.")
    # Calculate A^3
    A_cubed = np.linalg.matrix_power(A, 3)
    
    # Calculate A* (the conjugate transpose, or adjoint)
    A_star = A.conj().T

    print("A^3 =")
    print(A_cubed)
    print("A* =")
    print(A_star)

    # Check for equality using numpy.allclose for floating point precision
    is_satisfied = np.allclose(A_cubed, A_star)
    
    if is_satisfied:
        print("\nThe condition A^3 = A* is satisfied.")
    else:
        print("\nError: The condition A^3 = A* is NOT satisfied.")
        return

    print("\nStep 4: Find the set S of non-real eigenvalues for our matrix A and get its size.")
    
    # Get all eigenvalues from the matrix A
    eigenvalues_of_A = np.linalg.eigvals(A)
    
    # Filter for the ones that are not real (have a non-zero imaginary part)
    S = {val for val in eigenvalues_of_A if val.imag != 0}

    # The final equation is |S| = count
    # Outputting the numbers in the final equation as requested
    final_equation_lhs = "|S|"
    final_equation_rhs = len(S)
    
    print("The set of non-real eigenvalues S is:", S)
    print(f"The final equation is: {final_equation_lhs} = {final_equation_rhs}")

solve_eigenvalue_problem()