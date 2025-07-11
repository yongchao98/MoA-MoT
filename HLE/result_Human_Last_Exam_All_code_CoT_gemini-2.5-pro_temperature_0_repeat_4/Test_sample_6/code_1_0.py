import numpy as np

def count_kk_modes():
    """
    Calculates the number of Kaluza-Klein spin-2 mode eigenvalues below a threshold.

    This function solves the eigenvalue problem for KK modes in a warped compactification
    using the finite difference method.
    """
    # 1. Set up the numerical parameters
    N = 500  # Number of discretization points, chosen for good accuracy
    domain_length = 2 * np.pi
    dx = domain_length / N
    threshold = 14

    # 2. Define the grid and the potential term A'(x)
    # The grid runs from 0 to 2*pi, excluding the endpoint to handle periodicity.
    x_grid = np.linspace(0, domain_length, N, endpoint=False)
    
    # A(x) = sin(x) + 4*cos(x)
    # A'(x) = cos(x) - 4*sin(x)
    A_prime = np.cos(x_grid) - 4 * np.sin(x_grid)

    # 3. Construct the finite difference matrix for the operator
    # The ODE is: psi'' + 3*A'*psi' + m^2*psi = 0
    # Or: -psi'' - 3*A'*psi' = m^2*psi
    # We discretize the operator on the left-hand side.
    # Using central differences:
    # psi'(x_i)  ~ (psi_{i+1} - psi_{i-1}) / (2*dx)
    # psi''(x_i) ~ (psi_{i+1} - 2*psi_i + psi_{i-1}) / dx^2
    #
    # The matrix equation is M * psi = m^2 * psi.
    # M[i,i] * psi_i + M[i,i+1] * psi_{i+1} + M[i,i-1] * psi_{i-1} = m^2 * psi_i
    #
    # -( (psi_{i+1} - 2*psi_i + psi_{i-1})/dx^2 ) - 3*A'_i*( (psi_{i+1} - psi_{i-1})/(2*dx) ) = m^2 * psi_i
    #
    # Grouping terms by psi_i, psi_{i+1}, psi_{i-1}:
    # (2/dx^2)*psi_i + (-1/dx^2 - 3*A'_i/(2*dx))*psi_{i+1} + (-1/dx^2 + 3*A'_i/(2*dx))*psi_{i-1} = m^2*psi_i
    
    M = np.zeros((N, N))

    # Coefficients from the discretization
    coeff_diag = 2 / dx**2
    coeff_plus1 = -1 / dx**2 - 3 * A_prime / (2 * dx)
    coeff_minus1 = -1 / dx**2 + 3 * A_prime / (2 * dx)

    # Populate the matrix M
    for i in range(N):
        M[i, i] = coeff_diag
        # Periodic boundary conditions are handled by the modulo operator %
        M[i, (i + 1) % N] = coeff_plus1[i]
        M[i, (i - 1 + N) % N] = coeff_minus1[i]

    # 4. Calculate the eigenvalues of the matrix
    eigenvalues = np.linalg.eigvals(M)

    # 5. Process the eigenvalues
    # The operator is self-adjoint, so eigenvalues must be real.
    # We take the real part to discard small imaginary noise from numerical computation.
    m_squared = np.real(eigenvalues)
    
    # Sort the eigenvalues for easier inspection
    m_squared.sort()

    # 6. Count the eigenvalues below the threshold
    eigenvalues_below_threshold = m_squared[m_squared < threshold]
    count = len(eigenvalues_below_threshold)

    # 7. Print the results as requested
    print("The equation for the Kaluza-Klein mass-squared m^2 is the eigenvalue equation:")
    print("psi''(x) + 3*(cos(x) - 4*sin(x))*psi'(x) + m^2*psi(x) = 0")
    print("with periodic boundary conditions on x in [0, 2*pi].\n")
    
    print(f"The eigenvalues m^2 found below the threshold of {threshold} are:")
    for eig_val in eigenvalues_below_threshold:
        # The prompt asks to "output each number in the final equation".
        # We interpret this as showing the values of m^2 that satisfy the condition.
        print(f"m^2 = {eig_val:.4f}")

    print(f"\nThe total number of eigenvalues below {threshold} is {count}.")
    
    # Final answer in the required format
    print(f"\n<<<{count}>>>")

if __name__ == '__main__':
    count_kk_modes()