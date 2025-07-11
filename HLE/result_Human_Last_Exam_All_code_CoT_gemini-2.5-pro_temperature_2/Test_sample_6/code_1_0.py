import numpy as np
from scipy.linalg import eigh

def count_kk_modes():
    """
    Calculates the number of Kaluza-Klein spin-2 eigenvalues below a threshold.
    """
    # 1. Setup the problem
    # The number of points for discretization. A larger number gives higher accuracy.
    N = 2000
    # The internal dimension is a circle with x in [0, 2*pi]
    x_max = 2.0 * np.pi
    
    # Create the discretized grid.
    # endpoint=False ensures periodicity (x[N] is not included, as it's identified with x[0]).
    x, h = np.linspace(0, x_max, N, endpoint=False, retstep=True)
    
    # 2. Define the warp factor and its derivatives
    A = np.sin(x) + 4.0 * np.cos(x)
    A_prime = np.cos(x) - 4.0 * np.sin(x)
    A_double_prime = -np.sin(x) - 4.0 * np.cos(x)
    
    # 3. Define the potential V(x) and weight function W(x)
    # This comes from transforming the original ODE to a Schrödinger-like equation.
    V = 1.5 * A_double_prime + (9.0 / 4.0) * A_prime**2
    W = np.exp(2.0 * A)
    
    # 4. Construct the matrices for the generalized eigenvalue problem M*f = lambda*B*f
    
    # 'M' represents the operator -d^2/dx^2 + V(x) with periodic boundary conditions.
    M = np.zeros((N, N))
    main_diag = 2.0 / h**2 + V
    off_diag = -1.0 / h**2
    
    np.fill_diagonal(M, main_diag)
    M += np.diag(np.full(N - 1, off_diag), k=1)
    M += np.diag(np.full(N - 1, off_diag), k=-1)
    
    # Corner elements enforce periodic boundary conditions
    M[0, N - 1] = off_diag
    M[N - 1, 0] = off_diag
    
    # 'B' is the diagonal matrix for the weight function W(x)
    B = np.diag(W)
    
    # 5. Solve the generalized eigenvalue problem
    # eigh is suitable for real symmetric matrices (M) and positive definite matrices (B),
    # and it conveniently returns eigenvalues sorted in ascending order.
    # The eigenvalues are the masses squared, m^2.
    mass_squared, _ = eigh(M, b=B)
    
    # 6. Count the eigenvalues below the threshold and print the results.
    threshold = 14
    count = 0
    
    print(f"Searching for eigenvalues m^2 such that m^2 < {threshold}:\n")
    
    # Iterate through the calculated eigenvalues and check against the threshold.
    # The question requires printing the comparison for each found eigenvalue.
    for val in mass_squared:
        if val < threshold:
            print(f"Found eigenvalue m^2 = {val:.5f}. Since {val:.5f} < {threshold}, we count it.")
            count += 1
        else:
            # Since eigenvalues are sorted, we can stop once we are above the threshold.
            break
            
    print(f"\nThe total number of eigenvalues below {threshold} is {count}.")
    
    # Output the final answer in the requested format
    print(f"\n<<<{count}>>>")

# Execute the function to get the answer
count_kk_modes()