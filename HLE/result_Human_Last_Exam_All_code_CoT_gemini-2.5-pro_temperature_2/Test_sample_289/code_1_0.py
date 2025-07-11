import cmath
import numpy as np

def solve_eigenvalue_problem():
    """
    This function solves for the largest possible size of a set of non-real
    eigenvalues for a matrix A satisfying A^3 = A^*.
    """
    print("The problem is to find the maximum number of distinct non-real eigenvalues")
    print("for a complex matrix A that satisfies the equation A^3 = A*.")
    print("\n--- Step 1: Derive the equation for eigenvalues ---")
    print("A matrix A satisfying A^3 = A* must be normal (AA* = A*A).")
    print("This allows us to reduce the matrix equation to a scalar equation for its eigenvalues lambda.")
    print("The resulting eigenvalue equation is: lambda^3 = conjugate(lambda).")
    print("Here, the numbers in the equation are the power, 3, and the implicit power on the conjugate, 1.")
    
    print("\n--- Step 2: Solve the eigenvalue equation ---")
    print("We solve lambda^3 = conjugate(lambda) by writing lambda in polar form: lambda = r * e^(i*theta).")
    print("This yields two conditions:")
    print("  1. For the magnitude r: r^3 = r  => r=0 or r=1.")
    print("  2. For the angle theta (when r=1): e^(i*3*theta) = e^(-i*theta) => 4*theta = 2*k*pi.")
    print("     This gives theta = k*pi/2 for k=0, 1, 2, 3.")

    possible_eigenvalues = set()

    # Case r=0 gives lambda = 0
    lambda_0 = complex(0, 0)
    possible_eigenvalues.add(lambda_0)
    print(f"\nSolution for r=0 gives the real eigenvalue: {lambda_0}")

    # Case r=1 gives eigenvalues on the unit circle
    print("Solutions for r=1 give the following eigenvalues for k=0,1,2,3:")
    k_values = [0, 1, 2, 3]
    for k in k_values:
        theta = k * cmath.pi / 2
        # Calculate eigenvalue: e^(i*theta)
        eigenvalue = cmath.exp(complex(0, theta))
        # Round to handle floating point inaccuracies for clean output
        eigenvalue = complex(round(eigenvalue.real, 5), round(eigenvalue.imag, 5))
        possible_eigenvalues.add(eigenvalue)
        print(f"  k={k}: theta={k}*pi/2 => lambda = {eigenvalue}")
        
    print("\n--- Step 3: Identify non-real eigenvalues ---")
    print(f"The complete set of possible eigenvalues is: {possible_eigenvalues}")
    
    non_real_eigenvalues = {val for val in possible_eigenvalues if val.imag != 0}
    
    print(f"The subset S of possible non-real eigenvalues is: {non_real_eigenvalues}")
    
    print("\n--- Step 4: Determine the largest size of S ---")
    max_size = len(non_real_eigenvalues)
    print("The largest possible size of the set of non-real eigenvalues is the number of")
    print("distinct non-real solutions we found.")
    print(f"Result: The largest size |S| is {max_size}.")

    print("\n--- Step 5: Verification with an example ---")
    if max_size > 0:
        print(f"We construct a matrix A with these {max_size} non-real eigenvalues to show this size is achievable.")
        
        diag_elements = list(non_real_eigenvalues)
        A = np.diag(np.array(diag_elements, dtype=complex))

        # Calculate A^3 and A*
        A_cubed = np.linalg.matrix_power(A, 3)
        A_star = A.conj().T

        print("Let A be the diagonal matrix of the non-real eigenvalues:")
        print("A =\n", A)
        print("\nWe check if A^3 = A^* for this matrix.")
        print("A^3 =\n", A_cubed)
        print("A^* (adjoint of A) =\n", A_star)

        if np.allclose(A_cubed, A_star):
            print("\nThe condition A^3 = A^* holds true.")
            print(f"This matrix has a set of non-real eigenvalues of size {max_size}, confirming our result.")
        else:
            print("\nVerification failed. An error occurred in the logic or calculations.")

if __name__ == '__main__':
    solve_eigenvalue_problem()
