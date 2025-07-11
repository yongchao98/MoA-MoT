import numpy as np

def solve_quantum_bootstrap():
    """
    This function implements the quantum bootstrap method to find the ground state energy (E)
    and the expectation value of x^2 for the potential V(x) = x^2 + x^4, with K=7.
    """
    # Define the search grid for Energy (E) and <x^2> (x2).
    # The ranges are chosen based on physical estimates for similar potentials.
    # The step size of 0.001 is fine enough for the required 3-digit precision.
    e_range = np.arange(1.0, 1.2, 0.001)
    x2_range = np.arange(0.3, 0.5, 0.001)

    min_E_found = -1
    min_x2_found = -1

    # Loop over E first, from low to high, to find the minimum energy.
    for E in e_range:
        # For each E, search for a valid <x^2>.
        for x2 in x2_range:
            # y is an array to store moments y_k = <x^{2k}>. We need up to y_7=<x^14>.
            y = np.zeros(8)
            y[0] = 1.0  # y_0 = <x^0> = 1
            y[1] = x2   # y_1 = <x^2> is a parameter

            # Calculate y_2 = <x^4> as a base case for the recursion.
            # This comes from the recurrence relation with k=1.
            y[2] = (-8 * y[1] + 4 * E * y[0]) / 12.0

            # Use the recurrence relation to calculate higher moments y_3 to y_7.
            # y_{k+1} = ( -8k*y_k + (8k-4)E*y_{k-1} + (2k-1)(2k-2)(2k-3)y_{k-2} ) / (8k+4)
            for k in range(2, 7): # Loop for k from 2 to 6
                numerator = -8*k*y[k] + (8*k-4)*E*y[k-1] + (2*k-1)*(2*k-2)*(2*k-3)*y[k-2]
                denominator = 8*k + 4
                y[k+1] = numerator / denominator

            # Construct the two moment matrices.
            # M_even_{ij} = y_{i+j} for i,j = 0,1,2,3.
            # M_odd_{ij}  = y_{i+j+1} for i,j = 0,1,2,3.
            M_even = np.array([[y[i+j] for j in range(4)] for i in range(4)])
            M_odd = np.array([[y[i+j+1] for j in range(4)] for i in range(4)])

            # Check if both matrices are positive semidefinite by checking if their
            # minimum eigenvalue is non-negative.
            try:
                eig_even = np.linalg.eigvalsh(M_even)
                eig_odd = np.linalg.eigvalsh(M_odd)
            except np.linalg.LinAlgError:
                # If matrix is not valid (e.g., contains NaNs), skip.
                continue

            # A small negative tolerance handles floating point inaccuracies.
            if np.min(eig_even) >= -1e-9 and np.min(eig_odd) >= -1e-9:
                min_E_found = E
                min_x2_found = x2
                # We found the first valid <x^2> for this E. Break the inner loop.
                break
        
        if min_E_found != -1:
            # Since we iterate E from low to high, the first E that has a valid
            # <x^2> is the minimum E. We can stop the search.
            break

    # Print the final result.
    if min_E_found != -1:
        print("The minimal value for E and the corresponding value for <x^2> are found:")
        print(f"E = {min_E_found:.3f}")
        print(f"<x^2> = {min_x2_found:.3f}")
    else:
        print("No solution was found in the specified search range.")

# Run the simulation.
solve_quantum_bootstrap()