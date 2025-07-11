import numpy as np

def find_minimal_values_for_anharmonic_oscillator():
    """
    Solves for the ground state energy (E) and expectation value <x^2> of the
    anharmonic oscillator with potential V(x) = x^2 + x^4 using the
    quantum bootstrap method with K=7.
    """

    # Step 1 & 2: Define a function to calculate moments <x^n> using the recursion relation
    def calculate_moments(E, m2, K):
        """
        Calculates even moments <x^{2n}> up to n=K using the recursion relation.
        E: Energy
        m2: Expectation value <x^2>
        K: The truncation order of the operator O = sum(c_i * x^i)
        """
        # We need moments up to <x^(2K)>. For K=7, this is <x^14>.
        # Our moment array `m` will store m[i] = <x^i>.
        max_moment_idx = 2 * K
        m = np.zeros(max_moment_idx + 1)

        # Step 4: Initial conditions from the problem description
        m[0] = 1.0  # <x^0> = 1
        m[2] = m2   # <x^2> is a parameter we search over

        # The recursion relation derived from the bootstrap conditions for V(x) = x^2 + x^4 is:
        # (4t + 8) <x^{t+3}> = 4tE <x^{t-1}> + t(t-1)(t-2) <x^{t-3}> - (4t + 4) <x^{t+1}>
        # We use odd values of t to connect the even-powered moments.
        for t in range(1, max_moment_idx, 2):
            if t + 3 > max_moment_idx:
                break

            # Handle the <x^{t-3}> term for t=1, where its coefficient t(t-1)(t-2) is 0.
            m_t_minus_3 = 0 if t < 3 else m[t - 3]
                
            numerator = (4 * t * E * m[t - 1] +
                         t * (t - 1) * (t - 2) * m_t_minus_3 -
                         (4 * t + 4) * m[t + 1])
            denominator = 4 * t + 8
            
            m[t + 3] = numerator / denominator

        return m

    # Step 3 & 4: Check the positive semidefinite constraint on the moment matrix
    def check_positive_semidefinite(E, m2, K):
        """
        Calculates moments for a given E and <x^2>, constructs the moment matrices,
        and checks if they are positive semidefinite by examining their eigenvalues.
        """
        moments = calculate_moments(E, m2, K)

        # The moment matrix M_ij = <x^(i+j)> for i,j in {0..K} decouples
        # into an even and an odd block for a symmetric potential.
        
        # For K=7, the 8x8 matrix splits into two 4x4 matrices.
        # M_even_{ij} = <x^{2i+2j}> for i,j in {0,1,2,3}
        M_even = np.array([
            [moments[0], moments[2], moments[4], moments[6]],
            [moments[2], moments[4], moments[6], moments[8]],
            [moments[4], moments[6], moments[8], moments[10]],
            [moments[6], moments[8], moments[10], moments[12]]
        ])

        # M_odd_{ij} = <x^{2(i+j)+2}> for i,j in {0,1,2,3}
        M_odd = np.array([
            [moments[2], moments[4], moments[6], moments[8]],
            [moments[4], moments[6], moments[8], moments[10]],
            [moments[6], moments[8], moments[10], moments[12]],
            [moments[8], moments[10], moments[12], moments[14]]
        ])
        
        # A matrix is positive semidefinite if all its eigenvalues are non-negative.
        # We use a small tolerance for floating point inaccuracies.
        try:
            eig_even = np.linalg.eigvalsh(M_even)
            eig_odd = np.linalg.eigvalsh(M_odd)
        except np.linalg.LinAlgError:
            return False # Fails check if matrix is not valid for eigenvalue calculation

        return np.all(eig_even >= -1e-9) and np.all(eig_odd >= -1e-9)

    # Step 5: Search for the minimal E and corresponding <x^2>
    
    # Set the truncation order as specified in the problem
    K_val = 7
    
    # Based on known results, the ground state energy is ~1.392.
    # We search in a grid around this value with a step of 0.001 for 3-digit precision.
    E_values = np.arange(1.3, 1.5, 0.001)
    m2_values = np.arange(0.3, 0.5, 0.001)

    # We iterate E from low to high. The first E with a valid <x^2> is our minimum.
    for E in E_values:
        # For that E, we iterate <x^2> from low to high to find the minimum valid value.
        for m2 in m2_values:
            if check_positive_semidefinite(E, m2, K_val):
                # We have found the solution: the minimum E and the minimum
                # corresponding <x^2> on our search grid.
                print("Found the minimal values satisfying the bootstrap constraints:")
                print(f"Minimal E: {E:.3f}")
                print(f"Minimal <x^2>: {m2:.3f}")
                return

    print("No solution was found in the specified search range.")
    print("Consider expanding the search grid for E and <x^2>.")

if __name__ == '__main__':
    find_minimal_values_for_anharmonic_oscillator()