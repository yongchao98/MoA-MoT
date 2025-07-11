import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_kk_masses():
    """
    This script calculates the number of Kaluza-Klein graviton modes
    with mass-squared below 14 for the given 5D warped compactification.
    """

    # Define the warp factor and its derivative from the problem statement.
    # A(x) = sin(x) + 4*cos(x)
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    def A_prime(x):
        return np.cos(x) - 4 * np.sin(x)

    # Define the system of first-order ODEs.
    # We solve for two independent solutions simultaneously.
    # y = [psi_1, psi_1', psi_2, psi_2']
    def ode_system(x, y, m_squared):
        # The ODE is psi'' + 3*A'*psi' + m^2*exp(2A)*psi = 0
        psi1_pp = -3 * A_prime(x) * y[1] - m_squared * np.exp(2 * A(x)) * y[0]
        psi2_pp = -3 * A_prime(x) * y[3] - m_squared * np.exp(2 * A(x)) * y[2]
        return [y[1], psi1_pp, y[3], psi2_pp]

    # This function solves the ODE for a given m_squared and returns the
    # final state vector y(2*pi), which defines the monodromy matrix M.
    # M = [[y[0], y[2]], [y[1], y[3]]] at x=2*pi.
    def get_final_y(m_squared):
        y0 = [1.0, 0.0, 0.0, 1.0] # Initial conditions for two basis solutions
        x_span = [0, 2 * np.pi]
        
        # High precision is needed for stable results
        sol = solve_ivp(
            lambda x, y: ode_system(x, y, m_squared),
            x_span,
            y0,
            rtol=1e-9,
            atol=1e-9,
            method='Radau' # Good for stiff problems
        )
        return sol.y[:, -1]

    # The objective function for the root finder.
    # An eigenvalue exists when Tr(M) = y1(2pi) + y2'(2pi) = 2.
    def trace_objective_func(m_squared):
        if m_squared < 0: return 1e9 # Mass squared must be non-negative
        y_final = get_final_y(m_squared)
        trace = y_final[0] + y_final[3]
        return trace - 2

    # --- Main Calculation ---
    
    max_m_squared = 14.0
    eigenvalues = []

    # The m^2=0 mode (4D graviton) is always an eigenvalue.
    # It is non-degenerate for this warp factor.
    eigenvalues.append({'value': 0.0, 'degenerate': False})

    # Search for other eigenvalues by finding roots of trace_objective_func
    # We create a grid to find intervals where a root exists (sign change).
    grid_points = 1000
    m2_grid = np.linspace(0.01, max_m_squared, grid_points)
    
    f_prev = trace_objective_func(m2_grid[0])

    for i in range(1, len(m2_grid)):
        m2_curr = m2_grid[i]
        f_curr = trace_objective_func(m2_curr)
        
        # A sign change indicates a root is bracketed
        if f_prev * f_curr < 0:
            m2_prev = m2_grid[i-1]
            try:
                # Find the precise root
                root = brentq(trace_objective_func, m2_prev, m2_curr)
                
                # Check for degeneracy. This occurs if M=I.
                # A sufficient check is if the off-diagonal elements are zero.
                # We check y2(2pi).
                y_final_at_root = get_final_y(root)
                y2_at_2pi = y_final_at_root[2]
                
                # A very small value of y2(2pi) implies degeneracy.
                is_degenerate = np.isclose(y2_at_2pi, 0.0, atol=1e-7)
                
                eigenvalues.append({'value': root, 'degenerate': is_degenerate})
                
            except ValueError:
                # This may happen if function values are not finite,
                # but should be avoided by robust integration.
                pass

        f_prev = f_curr
        
    # --- Output the results ---
    
    total_count = 0
    print("Found eigenvalues m^2 below 14:")
    for ev in sorted(eigenvalues, key=lambda x: x['value']):
        total_count += 2 if ev['degenerate'] else 1
        degeneracy_str = "degenerate" if ev['degenerate'] else "non-degenerate"
        print(f"  m^2 = {ev['value']:.4f} ({degeneracy_str})")

    print("\nCounting these eigenvalues (with degeneracy):")
    equation_parts = []
    for ev in sorted(eigenvalues, key=lambda x: x['value']):
        if ev['degenerate']:
            equation_parts.append("2")
        else:
            equation_parts.append("1")
    
    final_equation = " + ".join(equation_parts)
    print(f"Total Count = {final_equation} = {total_count}")
    
    return total_count

# Run the calculation and print the final answer in the required format
final_count = solve_kk_masses()
print(f"<<<{final_count}>>>")
