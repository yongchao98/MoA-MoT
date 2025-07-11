import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_for_kk_modes():
    """
    This script finds the number of Kaluza-Klein eigenvalues for a given
    5D gravitational theory compactified on a circle, below a certain threshold.
    """
    eigenvalue_threshold = 14
    
    # 1. Define the warp factor and its derivative from the problem statement.
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    def A_prime(x):
        return np.cos(x) - 4 * np.sin(x)

    # 2. Define the ODE system y' = f(x, y) where y = [phi, phi'].
    # The ODE is phi'' + 3*A'(x)*phi' + m^2*exp(2*A(x))*phi = 0.
    def ode_system(x, y, m_squared):
        phi, phi_prime = y
        d_phi_dx = phi_prime
        d_phi_prime_dx = -3 * A_prime(x) * phi_prime - m_squared * np.exp(2 * A(x)) * phi
        return [d_phi_dx, d_phi_prime_dx]

    # Memoization cache for the trace function to avoid re-computation.
    memo_trace = {}

    # 3. Define a function to compute the trace of the monodromy matrix for a given m^2.
    def get_trace(m_squared):
        """
        Calculates the trace of the monodromy matrix by solving the ODE
        for two independent initial conditions.
        """
        if m_squared in memo_trace:
            return memo_trace[m_squared]

        x_span = [0, 2 * np.pi]
        
        # Solution for initial condition y(0) = [1, 0]
        sol1 = solve_ivp(ode_system, x_span, [1, 0], args=(m_squared,), rtol=1e-8, atol=1e-8)
        a, b = sol1.y[:, -1]
        
        # Solution for initial condition y(0) = [0, 1]
        sol2 = solve_ivp(ode_system, x_span, [0, 1], args=(m_squared,), rtol=1e-8, atol=1e-8)
        c, d = sol2.y[:, -1]
        
        trace = a + d
        memo_trace[m_squared] = trace
        return trace

    # The objective function whose roots give the eigenvalues m^2.
    # An eigenvalue exists when trace(M) = 2.
    def objective_func(m_squared):
        return get_trace(m_squared) - 2

    # 4. Count the eigenvalues below the threshold.
    found_positive_roots = []
    
    # Start with the non-degenerate zero-mode m^2 = 0.
    total_eigenvalues = 1
    print(f"Found eigenvalue m^2 = 0. This mode is non-degenerate. Running count: {total_eigenvalues}")

    # Scan for the positive, doubly-degenerate eigenvalues.
    # We search for roots of objective_func in small intervals.
    search_step = 0.5
    m_sq_low = 1e-6  # Start just above 0 to avoid the root we already counted.

    while m_sq_low < eigenvalue_threshold:
        m_sq_high = m_sq_low + search_step
        if m_sq_high > eigenvalue_threshold:
            m_sq_high = eigenvalue_threshold
        
        # A root exists in (low, high) if the function values at the endpoints have different signs.
        try:
            val_low = objective_func(m_sq_low)
            val_high = objective_func(m_sq_high)

            if np.sign(val_low) != np.sign(val_high):
                # Use Brent's method to find the root accurately.
                root = brentq(objective_func, m_sq_low, m_sq_high)
                
                # Add the root if it's within the threshold and not already found
                if root < eigenvalue_threshold and not np.isclose(root, found_positive_roots[-1:]):
                    found_positive_roots.append(root)
                    total_eigenvalues += 2 # Each positive root corresponds to a doubly-degenerate state.
                    print(f"Found eigenvalue m^2 = {root:.4f}. This mode is doubly degenerate. Running count: {total_eigenvalues}")
        except (ValueError, IndexError):
            # ValueError: brentq fails if signs are the same.
            # IndexError: handles the initial empty list case.
            pass
        
        # Move to the next search interval.
        m_sq_low = m_sq_high
        
        if m_sq_low >= eigenvalue_threshold:
            break

    # 5. Print the final result and the summation equation.
    print("\n-------------------------------------------------------------")
    print(f"Found {len(found_positive_roots)} doubly-degenerate positive eigenvalues below m^2 = {eigenvalue_threshold}.")
    
    equation_parts = ["1"] # Start with the m^2=0 mode
    for _ in found_positive_roots:
        equation_parts.append("2")
        
    equation = " + ".join(equation_parts)
    print(f"The total number of eigenvalues is given by the sum:")
    print(f"{equation} = {total_eigenvalues}")
    print("-------------------------------------------------------------")
    
    return total_eigenvalues

# Execute the function to solve the problem
solve_for_kk_modes()