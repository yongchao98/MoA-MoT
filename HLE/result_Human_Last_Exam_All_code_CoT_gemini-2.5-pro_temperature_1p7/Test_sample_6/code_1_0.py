import numpy as np
from scipy.linalg import eigh

def solve_graviton_masses():
    """
    This script calculates the number of Kaluza-Klein graviton mass eigenvalues
    below a certain threshold for the specified 5D warped compactification.
    It discretizes the corresponding Sturm-Liouville problem and solves it
    as a generalized matrix eigenvalue problem.
    """
    
    # Set up numerical parameters for the discretization.
    # A sufficiently large number of grid points is used for accuracy.
    N = 500  # Number of grid points
    L = 2 * np.pi  # Period of the compact dimension
    THRESHOLD = 14.0  # The upper bound for the eigenvalues m^2 to be counted.
    
    h = L / N  # Grid spacing
    # Create the grid for the coordinate x, from 0 to L-h.
    x = np.linspace(0, L, N, endpoint=False)

    # Define the warp factor A(x) and the derived functions p(x) and w(x)
    # for the Sturm-Liouville problem: -(p(x)psi'(x))' = m^2 * w(x) * psi(x).
    def A_func(t):
        return np.sin(t) + 4 * np.cos(t)

    def p_func(t):
        # p(x) = exp(3*A(x))
        return np.exp(3 * A_func(t))

    def w_func(t):
        # w(x) = exp(A(x))
        return np.exp(A_func(t))

    # Construct the matrices M and W for the generalized eigenvalue problem M*v = lambda*W*v.
    
    # W is a diagonal matrix containing the values of the weight function w(x) at each grid point.
    w_values = w_func(x)
    W = np.diag(w_values)

    # M is the symmetric matrix representing the discretized operator -(d/dx p(x) d/dx).
    # For better accuracy, we evaluate p(x) at the midpoints between grid points.
    x_midpoints = x + h / 2
    p_midpoints = p_func(x_midpoints)
    
    # Initialize the matrix M.
    M = np.zeros((N, N))
    
    # Fill the matrix M according to the finite difference scheme with periodic boundary conditions.
    for i in range(N):
        i_minus_1 = (i - 1 + N) % N # Index of the previous point with periodicity
        
        # Get the values of p at the midpoints surrounding x_i
        p_i_plus_half = p_midpoints[i]
        p_i_minus_half = p_midpoints[i_minus_1]
        
        # Main diagonal element
        M[i, i] = (p_i_plus_half + p_i_minus_half) / h**2
        
        # Off-diagonal elements, ensuring symmetry M[i,j] = M[j,i]
        i_plus_1 = (i + 1) % N # Index of the next point with periodicity
        M[i, i_plus_1] = -p_i_plus_half / h**2
        M[i, i_minus_1] = -p_i_minus_half / h**2

    # Solve the generalized eigenvalue problem M*v = lambda*W*v.
    # We use eigh because M is symmetric and W is symmetric positive-definite.
    # The eigenvalues (lambda) returned are the squared masses m^2.
    eigenvalues, _ = eigh(M, W)

    # Filter the eigenvalues to find those below the specified threshold.
    eigenvalues_below_threshold = eigenvalues[eigenvalues < THRESHOLD]
    count = len(eigenvalues_below_threshold)

    # Output the results as requested.
    # First, print each eigenvalue found below the threshold.
    print("The eigenvalues m^2 below 14 are:")
    for ev in eigenvalues_below_threshold:
        print(f"{ev:.4f}")

    # Then, print the final count.
    print(f"\nTotal number of eigenvalues below {THRESHOLD} is {count}.")
    
    # Finally, print the answer in the required format.
    print(f"<<<{count}>>>")

# Run the calculation and print the output.
solve_graviton_masses()