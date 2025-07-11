import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_kaluza_klein_masses():
    """
    This script finds the number of Kaluza-Klein eigenvalues for spin-2 modes
    below a specified threshold in a warped compactification.

    The eigenvalue equation is: (e^{2A(x)} * phi')' + m^2 * phi = 0
    with A(x) = sin(x) + 4*cos(x) and periodic boundary conditions on [0, 2*pi].
    """
    
    # Define the warp factor and its derivative from the problem statement.
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    def A_prime(x):
        return np.cos(x) - 4 * np.sin(x)

    # Define the system of first-order ODEs for the shooting method.
    # y[0] = phi, y[1] = phi', lambda_val = m^2
    def ode_system(x, y, lambda_val):
        phi, phi_prime = y
        # The ODE is phi'' = -2*A'(x)*phi' - lambda_val*exp(-2*A(x))*phi
        d_phi_dx = phi_prime
        d_phi_prime_dx = -2 * A_prime(x) * phi_prime - lambda_val * np.exp(-2 * A(x)) * phi
        return [d_phi_dx, d_phi_prime_dx]

    # Define the characteristic function whose roots give the eigenvalues m^2.
    # An eigenvalue exists when Tr(M) - 2 = 0, where M is the monodromy matrix.
    def char_func(lambda_val):
        if lambda_val < 1e-6: # Handle m^2=0 case separately to avoid numerical issues
            return 0.0

        # High precision is needed for the numerical integration.
        rtol, atol = 1e-9, 1e-12
        
        # Integrate with initial condition (phi(0), phi'(0)) = (1, 0)
        sol1 = solve_ivp(ode_system, [0, 2 * np.pi], [1.0, 0.0], args=(lambda_val,), rtol=rtol, atol=atol)
        phi1_end, phi1_prime_end = sol1.y[:, -1]

        # Integrate with initial condition (phi(0), phi'(0)) = (0, 1)
        sol2 = solve_ivp(ode_system, [0, 2 * np.pi], [0.0, 1.0], args=(lambda_val,), rtol=rtol, atol=atol)
        phi2_end, phi2_prime_end = sol2.y[:, -1]
        
        # The trace of the monodromy matrix is phi1(2pi) + phi2'(2pi)
        trace_m = phi1_end + phi2_prime_end
        return trace_m - 2

    # --- Main Calculation ---
    
    # The m^2=0 mode (4D graviton) is always a solution and is non-degenerate.
    eigenvalue_count = 1
    eigenvalues = [0.0]
    
    threshold = 14.0

    # Scan for positive eigenvalues by finding roots of char_func.
    # A fine grid is used to not miss any roots.
    grid_points = 1500
    lambda_grid = np.linspace(0.01, threshold, grid_points)
    char_values = np.array([char_func(l) for l in lambda_grid])
    
    # Find intervals where the function changes sign, indicating a root.
    for i in range(len(char_values) - 1):
        if np.sign(char_values[i]) != np.sign(char_values[i+1]):
            a, b = lambda_grid[i], lambda_grid[i+1]
            try:
                # Find the root accurately within the interval.
                root = brentq(char_func, a, b)
                if root < threshold:
                    eigenvalues.append(root)
            except (ValueError, RuntimeError):
                # This should not happen given the sign change check, but is here for safety.
                pass
    
    # Count the eigenvalues, remembering that positive mass states are doubly degenerate.
    final_count = 1  # Start with the non-degenerate zero mode
    
    # Build the explanatory output string
    positive_eigenvalues = sorted([e for e in eigenvalues if e > 0])
    
    equation_parts = ["1 (for m^2=0.00)"]
    for val in positive_eigenvalues:
        final_count += 2
        equation_parts.append(f"2 (for m^2={val:.2f})")

    equation_str = " + ".join(equation_parts)
    print("Found eigenvalues and calculated the total count:")
    print(f"Number of eigenvalues = {equation_str} = {final_count}")
    
    return final_count

# Execute the calculation and store the final answer
final_answer = solve_kaluza_klein_masses()
# The final answer is wrapped in <<<>>> as requested.
# The code itself prints the detailed breakdown.
# print(f"<<<{final_answer}>>>")