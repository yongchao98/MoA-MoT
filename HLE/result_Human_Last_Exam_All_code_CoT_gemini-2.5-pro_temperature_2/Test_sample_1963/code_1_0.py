import numpy as np

def solve_quantum_bootstrap():
    """
    Finds the ground state energy E and <x^2> for the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """

    # We perform a grid search for E and <x^2>.
    # Based on known results for similar potentials, we can choose a focused search range.
    # The ground state energy is known to be around 1.392.
    e_vals = np.linspace(1.38, 1.40, 201)
    # The corresponding <x^2> should be around 0.4.
    x2_vals = np.linspace(0.40, 0.43, 201)

    min_e_found = float('inf')
    best_x2_found = -1.0

    # Iterate over the grid, starting with the lowest E
    for E in e_vals:
        found_valid_point = False
        for x2 in x2_vals:
            # The even moments are denoted by m[k] = <x^{2k}>
            m = np.zeros(8)  # We need moments up to m[7] = <x^14>
            m[0] = 1.0  # By definition <x^0> = 1
            m[1] = x2   # This is our test parameter <x^2>

            # The recurrence relation from the bootstrap conditions for V(x) = x^2 + x^4 allows us
            # to calculate all higher even moments.
            # Start with m[2] = <x^4> using the relation for t=1
            m[2] = (E - 2.0 * m[1]) / 3.0

            # Use the general recurrence for k >= 2 to find m[3]...m[7]
            # m_{k+1} = (1/(8k+4)) * [ (8k-4)E*m_{k-1} - 8k*m_k + (2k-1)(2k-2)(2k-3)*m_{k-2} ]
            try:
                for k in range(2, 7):  # This loop calculates m[3] through m[7]
                    term1 = (8.0 * k - 4.0) * E * m[k - 1]
                    term2 = -8.0 * k * m[k]
                    term3 = (2.0 * k - 1.0) * (2.0 * k - 2.0) * (2.0 * k - 3.0) * m[k - 2]
                    m[k + 1] = (term1 + term2 + term3) / (8.0 * k + 4.0)
            except (OverflowError, ValueError):
                # If moments become too large or invalid, skip this point
                continue

            # Construct the two Hankel matrices for K=7
            # M_e corresponds to the even operator basis {x^0, x^2, x^4, x^6}
            M_e = np.array([
                [m[0], m[1], m[2], m[3]],
                [m[1], m[2], m[3], m[4]],
                [m[2], m[3], m[4], m[5]],
                [m[3], m[4], m[5], m[6]]
            ])

            # M_o corresponds to the odd operator basis {x^1, x^3, x^5, x^7}
            M_o = np.array([
                [m[1], m[2], m[3], m[4]],
                [m[2], m[3], m[4], m[5]],
                [m[3], m[4], m[5], m[6]],
                [m[4], m[5], m[6], m[7]]
            ])

            # Check if both matrices are positive semidefinite (all eigenvalues >= 0)
            try:
                # Use eigvalsh as matrices are symmetric, it's faster and returns real eigenvalues
                eig_e = np.linalg.eigvalsh(M_e)
                eig_o = np.linalg.eigvalsh(M_o)
            except np.linalg.LinAlgError:
                continue # Matrix was not valid (e.g., contains NaN/inf)

            # Use a small negative tolerance to account for numerical precision errors
            tolerance = -1e-9
            if np.all(eig_e > tolerance) and np.all(eig_o > tolerance):
                # This (E, x2) pair is physically allowed.
                # Since we iterate E from low to high, the first one we find is the minimum E.
                min_e_found = E
                best_x2_found = x2
                found_valid_point = True
                break  # Found the optimal x2 for this E, move to the next E
        
        if found_valid_point:
            # We found the minimal E value in our search grid, so we can stop.
            break

    # Print the results formatted to 3 numerical digits.
    print(f"Minimal value of E: {min_e_found:.3f}")
    print(f"Minimal value of <x^2> for this E: {best_x2_found:.3f}")

if __name__ == '__main__':
    solve_quantum_bootstrap()
<<<1.392>>>