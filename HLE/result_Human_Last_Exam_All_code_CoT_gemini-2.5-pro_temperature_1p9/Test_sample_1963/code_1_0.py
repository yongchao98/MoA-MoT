import numpy as np

def calculate_moments(E, M1, max_n):
    """
    Calculates the moments M_n = <x^{2n}> using the recurrence relation.
    
    Args:
        E (float): The energy value.
        M1 (float): The value of <x^2>.
        max_n (int): The maximum index 'n' for M_n to calculate.

    Returns:
        numpy.ndarray: An array of moments [M_0, M_1, ..., max_n], or None on failure.
    """
    if M1 <= 0:
        return None
    
    moments = np.zeros(max_n + 1, dtype=np.float64)
    moments[0] = 1.0  # M_0 = <x^0> = 1
    moments[1] = M1

    # Use the recurrence M_{n+2} = f(M_{n+1}, M_n, M_{n-1}) derived from the bootstrap conditions.
    # We solve for <x^{t+3}>: <x^{t+3}> = (1/(4t+8)) * [4tE<x^{t-1}> + t(t-1)(t-2)<x^{t-3}> - (4t+4)<x^{t+1}>]
    # Set t = 2n+1 to get recurrence for M_n = <x^{2n}>
    # M_{n+2} = (1/(8n+12)) * [(8n+4)E*M_n + (2n+1)(2n)(2n-1)*M_{n-1} - (8n+8)*M_{n+1}]
    
    try:
        # Calculate M_2 (case n=0)
        # 12<x^4> = 4E - 8<x^2> => M_2 = E/3 - (2/3)*M_1
        moments[2] = E / 3.0 - (2.0 / 3.0) * M1

        # Calculate M_3, M_4, ... , M_{max_n} for n >= 1
        for n in range(1, max_n - 1):
            M_n_minus_1 = moments[n - 1]
            M_n = moments[n]
            M_n_plus_1 = moments[n + 1]

            term_M_n_minus_1_coeff = (2 * n + 1) * (2 * n) * (2 * n - 1)
            
            numerator = (8 * n + 4) * E * M_n + term_M_n_minus_1_coeff * M_n_minus_1 - (8 * n + 8) * M_n_plus_1
            denominator = 8 * n + 12
            
            if abs(denominator) < 1e-12: return None
            
            moments[n + 2] = numerator / denominator

    except (OverflowError, ValueError):
        return None # Return None if recurrence becomes numerically unstable
        
    return moments

def is_psd(matrix):
    """
    Checks if a matrix is positive semidefinite by examining its eigenvalues.
    A small tolerance is used for numerical stability.
    """
    if np.any(np.isnan(matrix)) or np.any(np.isinf(matrix)):
        return False
    # Use eigvalsh for symmetric matrices for better performance and numerical accuracy.
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.min(eigenvalues) >= -1e-9

def solve_bootstrap():
    """
    Performs the grid search to find the minimal E and corresponding <x^2>.
    """
    # K=7 requires constructing matrices from {1, x, ..., x^7}.
    # The highest moment needed is for the H_odd matrix, which requires M_7 = <x^14>.
    max_moment_n = 7 

    # Define the search grid. Based on physical arguments, E > 1 and <x^2> < 0.5.
    E_range = np.arange(1.2, 1.6, 0.001)
    M1_range = np.arange(0.3, 0.6, 0.001)

    min_E_found = float('inf')
    best_M1_found = -1

    for E in E_range:
        found_valid_M1 = False
        for M1 in M1_range:
            moments = calculate_moments(E, M1, max_moment_n)
            
            if moments is None:
                continue

            # Construct the 4x4 Hankel matrix for even powers
            H_even = np.array([
                [moments[0], moments[1], moments[2], moments[3]],
                [moments[1], moments[2], moments[3], moments[4]],
                [moments[2], moments[3], moments[4], moments[5]],
                [moments[3], moments[4], moments[5], moments[6]],
            ])

            # Construct the 4x4 Hankel matrix for odd powers
            H_odd = np.array([
                [moments[1], moments[2], moments[3], moments[4]],
                [moments[2], moments[3], moments[4], moments[5]],
                [moments[3], moments[4], moments[5], moments[6]],
                [moments[4], moments[5], moments[6], moments[7]],
            ])
            
            # If both matrices are positive semidefinite, we found a valid solution
            if is_psd(H_even) and is_psd(H_odd):
                min_E_found = E
                best_M1_found = M1
                found_valid_M1 = True
                break # Found a valid M1, move to the next E
        
        if found_valid_M1:
            # Since we iterate E from low to high, the first one found is the minimum
            break

    if best_M1_found != -1:
        print(f"Minimal E: {min_E_found:.3f}")
        print(f"Corresponding <x^2>: {best_M1_found:.3f}")
        
        # Recalculate and print the moments for the final answer
        final_moments = calculate_moments(min_E_found, best_M1_found, max_moment_n)
        print("\nCalculated moments M_n = <x^{2n}> for the solution:")
        for i, m in enumerate(final_moments):
            print(f"M_{i} = <x^{2*i}> = {m:.4f}")
    else:
        print("No solution found in the specified range.")

    # Return the final numerical answer in the specified format
    print(f"\n<<<Minimal E: {min_E_found:.3f}, Minimal <x^2>: {best_M1_found:.3f}>>>")

# Run the solver
solve_bootstrap()