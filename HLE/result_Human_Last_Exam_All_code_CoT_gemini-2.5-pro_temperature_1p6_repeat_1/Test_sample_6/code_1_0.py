import numpy as np
from scipy.integrate import solve_ivp

def solve_for_eigenvalues():
    """
    Calculates the number of Kaluza-Klein eigenvalues below 14 for the given
    5D gravitational theory.
    """
    # Define the warp factor and its derivative
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    def A_prime(x):
        return np.cos(x) - 4 * np.sin(x)

    # Define the potential in the Schrödinger equation
    # V(x) = (3/2)A'' + (9/4)(A')^2 = -(3/2)A + (9/4)(A')^2
    def V(x):
        return -1.5 * A(x) + 2.25 * (A_prime(x)**2)

    # Define the ODE system for the Schrödinger equation:
    # y[0] = phi, y[1] = phi_prime
    # phi'' = (V(x) - m^2) * phi
    def ode_system(x, y, m_squared):
        phi, phi_prime = y
        return [phi_prime, (V(x) - m_squared) * phi]

    # Set up the scan for m^2
    m_squared_limit = 14
    # A fine grid to avoid missing sign changes
    m_squared_vals = np.linspace(0.001, m_squared_limit, 2000)

    # --- Count Neumann eigenvalues (for even solutions) ---
    # Start with 1 to account for the m^2 = 0 eigenvalue (constant mode), which is an even solution.
    neumann_count = 1
    y0_neumann = [1.0, 0.0]  # IC: phi(0)=1, phi'(0)=0

    # Get the value at the start of the scan interval
    sol_N = solve_ivp(ode_system, [0, np.pi], y0_neumann, args=(m_squared_vals[0],), dense_output=True)
    # The value to check for sign changes is phi'(pi)
    prev_val_N = sol_N.sol(np.pi)[1]

    for m_squared in m_squared_vals[1:]:
        sol_N = solve_ivp(ode_system, [0, np.pi], y0_neumann, args=(m_squared,), dense_output=True)
        curr_val_N = sol_N.sol(np.pi)[1]
        # A sign change indicates an eigenvalue was crossed
        if np.sign(curr_val_N) != np.sign(prev_val_N):
            neumann_count += 1
        prev_val_N = curr_val_N

    # --- Count Dirichlet eigenvalues (for odd solutions) ---
    dirichlet_count = 0
    y0_dirichlet = [0.0, 1.0]  # IC: phi(0)=0, phi'(0)=1
    
    # Get the value at the start of the scan interval
    sol_D = solve_ivp(ode_system, [0, np.pi], y0_dirichlet, args=(m_squared_vals[0],), dense_output=True)
    # The value to check for sign changes is phi(pi)
    prev_val_D = sol_D.sol(np.pi)[0]

    for m_squared in m_squared_vals[1:]:
        sol_D = solve_ivp(ode_system, [0, np.pi], y0_dirichlet, args=(m_squared,), dense_output=True)
        curr_val_D = sol_D.sol(np.pi)[0]
        # A sign change indicates an eigenvalue was crossed
        if np.sign(curr_val_D) != np.sign(prev_val_D):
            dirichlet_count += 1
        prev_val_D = curr_val_D
    
    # The total count is the sum from both cases, accounting for degeneracies
    total_eigenvalues = neumann_count + dirichlet_count

    print(f"Number of eigenvalues (m^2 < 14) from even solutions (Neumann BCs): {neumann_count}")
    print(f"Number of eigenvalues (m^2 < 14) from odd solutions (Dirichlet BCs): {dirichlet_count}")
    print(f"Total number of eigenvalues below 14: {total_eigenvalues}")
    return total_eigenvalues

# Run the calculation and get the final answer.
final_answer = solve_for_eigenvalues()
print(f"<<<{final_answer}>>>")

solve_for_eigenvalues()