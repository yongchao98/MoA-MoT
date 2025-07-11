import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_kk_masses():
    """
    This script calculates the number of spin-2 Kaluza-Klein eigenvalues below 14
    for the given 5D gravitational theory.
    """

    # Define the warp factor A(x) and its derivative A'(x)
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    def A_prime(x):
        return np.cos(x) - 4 * np.sin(x)

    # Define the ODE system for y = [psi, psi']
    def ode_system(x, y, m2):
        psi, psi_prime = y
        # The ODE is: psi'' + 3*A'*psi' - m^2*exp(-2A)*psi = 0
        d_psi_prime_dx = -3 * A_prime(x) * psi_prime + m2 * np.exp(-2 * A(x)) * psi
        return [psi_prime, d_psi_prime_dx]

    # Function to compute the transfer matrix M for a given m^2
    def get_transfer_matrix(m2):
        x_span = [0, 2 * np.pi]
        
        # First solution with initial condition y(0) = [1, 0]
        sol1 = solve_ivp(ode_system, x_span, [1.0, 0.0], args=(m2,), rtol=1e-9, atol=1e-9)
        y11, y21 = sol1.y[:, -1]
        
        # Second solution with initial condition y(0) = [0, 1]
        sol2 = solve_ivp(ode_system, x_span, [0.0, 1.0], args=(m2,), rtol=1e-9, atol=1e-9)
        y12, y22 = sol2.y[:, -1]
        
        return np.array([[y11, y12], [y21, y22]])

    # The condition for an eigenvalue is Tr(M(m^2)) = 2.
    # We look for roots of the function f(m^2) = Tr(M(m^2)) - 2.
    def f(m2):
        if m2 < 0:
            return np.inf
        M = get_transfer_matrix(m2)
        return np.trace(M) - 2

    # Count the eigenvalues below 14.
    # Start with the massless mode m^2 = 0, which is non-degenerate.
    count = 1
    eigenvalues = [0.0]
    
    print("Searching for Kaluza-Klein eigenvalues m^2 < 14...")
    print(f"Found eigenvalue m^2 = {eigenvalues[0]:.4f} (non-degenerate, count = 1)")

    # Scan for massive (m^2 > 0) degenerate modes.
    # These are doubly degenerate, so each root adds 2 to the count.
    m2_step = 0.2
    m2_range = np.arange(m2_step, 14.0, m2_step)
    
    f_prev = f(1e-6) # Start with a value very close to 0
    
    found_roots = []
    
    for m2_curr in m2_range:
        f_curr = f(m2_curr)
        # Check for a sign change, which indicates a root in the interval
        if f_prev * f_curr < 0:
            try:
                root = brentq(f, m2_curr - m2_step, m2_curr)
                if root < 14.0:
                    found_roots.append(root)
                    count += 2
                    print(f"Found eigenvalue m^2 = {root:.4f} (degenerate, count = 2)")
            except ValueError:
                # brentq fails if f(a) and f(b) don't have opposite signs,
                # but we already checked for that. This is a safeguard.
                pass
        f_prev = f_curr

    # Final summary
    print("\n--- Calculation Summary ---")
    print("The total number of eigenvalues is the sum of counts for each mode found.")
    
    equation_parts = ["1"]
    for _ in found_roots:
        equation_parts.append("2")
    
    equation_str = " + ".join(equation_parts)
    print(f"Final count = {equation_str} = {count}")
    
    print("\nTotal number of eigenvalues below 14:")
    print(count)
    
    # Final answer in the required format
    print(f"<<<{count}>>>")

solve_kk_masses()