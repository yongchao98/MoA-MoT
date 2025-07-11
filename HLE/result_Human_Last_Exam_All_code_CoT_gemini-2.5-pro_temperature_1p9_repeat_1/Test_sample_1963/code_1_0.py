import numpy as np

def solve_bootstrap():
    """
    Solves for the ground state energy E and <x^2> of the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """

    # Define the search grid
    # Literature values suggest E_0 is around 1.39 and <x^2> is less than the QHO value of 0.5.
    E_values = np.arange(1.3, 1.5, 0.001)
    x2_values = np.arange(0.3, 0.5, 0.001)

    min_E_found = None
    min_x2_found = None
    final_moments = None

    # Grid search for the minimal allowed E and corresponding <x^2>
    for E in E_values:
        for x2 in x2_values:
            # Step 1: Calculate moments using the recursion relation
            moments = {0: 1.0, 2: x2}
            # Recursion for even moments m_2k is derived from the prompt's formula:
            # (8n+12)m_{2n+4} = -(8n+8)m_{2n+2} + 4(2n+1)E*m_{2n} + (2n+1)(2n)(2n-1)m_{2n-2}
            # We need moments up to m_14, which corresponds to n=5.
            try:
                for n in range(6):  # n from 0 to 5
                    m_2n = moments[2 * n]
                    m_2n_plus_2 = moments[2 * (n + 1)]
                    m_2n_minus_2 = moments.get(2 * (n - 1), 0.0)

                    numerator = (
                        -(8 * n + 8) * m_2n_plus_2
                        + 4 * (2 * n + 1) * E * m_2n
                        + (2 * n + 1) * (2 * n) * (2 * n - 1) * m_2n_minus_2
                    )
                    denominator = 8 * n + 12
                    
                    if denominator == 0:
                        raise ValueError("Division by zero in recursion.")

                    moments[2 * (n + 2)] = numerator / denominator

                # Step 2: Construct the moment matrices
                # M_ij = <x^(i+j)>. For K=7, i,j in {0..7}.
                # The matrix is block-diagonal. M_even for i,j even; M_odd for i,j odd.
                
                # M_even (indices 0, 2, 4, 6)
                M_even = np.zeros((4, 4))
                for i in range(4):
                    for j in range(4):
                        M_even[i, j] = moments[2*i + 2*j]

                # M_odd (indices 1, 3, 5, 7)
                M_odd = np.zeros((4, 4))
                for i in range(4):
                    for j in range(4):
                        # The moments are <x^( (2i+1) + (2j+1) )> = <x^(2i+2j+2)>
                        M_odd[i, j] = moments[2*i + 2*j + 2]

                # Step 3: Check if matrices are positive semidefinite
                eig_even = np.linalg.eigvalsh(M_even)
                eig_odd = np.linalg.eigvalsh(M_odd)

                if np.min(eig_even) >= -1e-9 and np.min(eig_odd) >= -1e-9:
                    min_E_found = E
                    min_x2_found = x2
                    final_moments = moments
                    break
            except (KeyError, np.linalg.LinAlgError):
                # Move to the next point if moments can't be computed or matrix is invalid
                continue
        
        if min_E_found is not None:
            break

    # Step 4: Output the results
    if min_E_found is not None:
        print(f"Minimal ground state energy E found: {min_E_found:.3f}")
        print(f"Corresponding expectation value <x^2> found: {min_x2_found:.3f}\n")

        print("Verifying the recursion relation for n=1 (t=3):")
        m0 = final_moments[0]
        m2 = final_moments[2]
        m4 = final_moments[4]
        m6 = final_moments[6]
        
        # The relation is: -20*<x^6> - 16*<x^4> + 12*E*<x^2> + 6*<x^0> = 0
        lhs = -20 * m6 - 16 * m4 + 12 * min_E_found * m2 + 6 * m0
        
        print("Equation: -20 * <x^6> - 16 * <x^4> + 12 * E * <x^2> + 6 * <x^0> = 0")
        print(f"Numerical check: -20 * {m6:.3f} - 16 * {m4:.3f} + 12 * {min_E_found:.3f} * {m2:.3f} + 6 * {m0:.3f} = {lhs:.4f}")
        
        # Final answer format
        # The prompt is ambiguous if it asks for one or two numbers.
        # Assuming it wants both minimal E and corresponding <x^2>.
        print(f"\n<<<E={min_E_found:.3f}, <x^2>={min_x2_found:.3f}>>>")
    else:
        print("Solution not found in the specified search range.")

solve_bootstrap()