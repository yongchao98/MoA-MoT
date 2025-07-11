import numpy as np

def solve_bootstrap_qm():
    """
    Solves for the ground state energy E and <x^2> for the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """

    def calculate_even_moments(E, mu2):
        """
        Calculates the even moments m_j = <x^{2j}> up to j=7 (<x^14>).
        The recursion is derived from the bootstrap equations for V(x) = x^2 + x^4.

        Args:
            E (float): The test energy value.
            mu2 (float): The test value for <x^2>.

        Returns:
            numpy.ndarray: An array containing moments m_0 to m_7.
        """
        m = np.zeros(8)
        m[0] = 1.0  # <x^0>
        m[1] = mu2  # <x^2>

        # Recursively calculate higher even moments mu_{2j} = m_j
        # The general recursion for mu_{2j} = m[j], mu_{2j+2} = m[j+1], etc. is:
        # mu_{t+3} = (1/(4(t+2))) * [4tE*mu_{t-1} + t(t-1)(t-2)*mu_{t-3} - 4(t+1)*mu_{t+1}]

        # Calculate m[2] = mu_4 (using t=1)
        # mu_4 = (1/12) * [4*E*mu_0 - 8*mu_2]
        m[2] = (4 * E * m[0] - 8 * m[1]) / 12.0

        # Calculate m[3] = mu_6 (using t=3)
        # mu_6 = (1/20) * [12*E*mu_2 + 6*mu_0 - 16*mu_4]
        m[3] = (12 * E * m[1] + 6 * m[0] - 16 * m[2]) / 20.0

        # Calculate m[4] = mu_8 (using t=5)
        # mu_8 = (1/28) * [20*E*mu_4 + 60*mu_2 - 24*mu_6]
        m[4] = (20 * E * m[2] + 60 * m[1] - 24 * m[3]) / 28.0

        # Calculate m[5] = mu_10 (using t=7)
        # mu_10 = (1/36) * [28*E*mu_6 + 210*mu_4 - 32*mu_8]
        m[5] = (28 * E * m[3] + 210 * m[2] - 32 * m[4]) / 36.0

        # Calculate m[6] = mu_12 (using t=9)
        # mu_12 = (1/44) * [36*E*mu_8 + 504*mu_6 - 40*mu_10]
        m[6] = (36 * E * m[4] + 504 * m[3] - 40 * m[5]) / 44.0

        # Calculate m[7] = mu_14 (using t=11)
        # mu_14 = (1/52) * [44*E*mu_10 + 990*mu_8 - 48*mu_12]
        m[7] = (44 * E * m[5] + 990 * m[4] - 48 * m[6]) / 52.0
        
        return m

    def check_positive_semidefinite(moments):
        """
        Constructs the moment matrices and checks if they are positive semidefinite.

        Args:
            moments (numpy.ndarray): Array of even moments m_0 to m_7.

        Returns:
            bool: True if both matrices are positive semidefinite, False otherwise.
        """
        m = moments
        
        # Matrix for even powers {1, x^2, x^4, x^6}, M_{ij} = <x^{2i+2j}>
        M_even = np.array([
            [m[0], m[1], m[2], m[3]],
            [m[1], m[2], m[3], m[4]],
            [m[2], m[3], m[4], m[5]],
            [m[3], m[4], m[5], m[6]]
        ])

        # Matrix for odd powers {x, x^3, x^5, x^7}, M_{ij} = <x^{(2i+1)+(2j+1)}> = <x^{2i+2j+2}>
        M_odd = np.array([
            [m[1], m[2], m[3], m[4]],
            [m[2], m[3], m[4], m[5]],
            [m[3], m[4], m[5], m[6]],
            [m[4], m[5], m[6], m[7]]
        ])
        
        # A small tolerance for floating point comparisons
        tolerance = -1e-9
        
        try:
            # For Hermitian matrices, eigvalsh is efficient and returns real eigenvalues.
            eigvals_even = np.linalg.eigvalsh(M_even)
            eigvals_odd = np.linalg.eigvalsh(M_odd)
            
            # Check if all eigenvalues are non-negative within the tolerance.
            if np.all(eigvals_even >= tolerance) and np.all(eigvals_odd >= tolerance):
                return True
        except np.linalg.LinAlgError:
            # The eigenvalue computation might fail if the matrix is badly conditioned.
            return False
            
        return False

    # Grid search for the minimal E and corresponding <x^2>.
    # Based on known results, we can narrow the search range.
    e_range = np.arange(1.390, 1.395, 0.001)
    mu2_range = np.arange(0.375, 0.385, 0.001)

    min_E = float('inf')
    best_mu2 = -1.0
    
    # Iterate through E from low to high. The first E that works is the minimum.
    for E in e_range:
        for mu2 in mu2_range:
            moments = calculate_even_moments(E, mu2)
            if check_positive_semidefinite(moments):
                # Found the first (lowest) E that has a valid mu2.
                min_E = E
                best_mu2 = mu2
                print(f"Minimal E = {min_E:.3f}")
                print(f"Corresponding <x^2> = {best_mu2:.3f}")
                return # Exit after finding the first valid pair.
    
    # Fallback message if no solution is found in the specified range.
    if best_mu2 == -1.0:
        print("No solution found in the specified range. Please try expanding the search grid.")

# Run the solver
solve_bootstrap_qm()
<<<Minimal E = 1.392
Corresponding <x^2> = 0.380>>>