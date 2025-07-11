import numpy as np

def calculate_susceptibility(c, J, beta, m0):
    """
    Calculates the magnetic susceptibility for an Ising model on a sparse random graph.

    Args:
        c (float): The connectivity (average degree) of the graph.
        J (float): The coupling constant.
        beta (float): The inverse temperature (1/T).
        m0 (float): The site magnetization.
    """
    # Check for physical validity of magnetization
    if not (0 <= m0 < 1):
        print("Error: Magnetization m0 must be in the range [0, 1).")
        return

    # 1. Calculate T_J = tanh(beta * J)
    T_J = np.tanh(beta * J)

    # 2. Calculate cavity magnetization m_c
    # artanh(m0) can lead to infinity if m0 is 1. We've checked for m0 < 1.
    artanh_m0 = np.arctanh(m0)
    m_c = np.tanh(((c - 1) / c) * artanh_m0)

    # 3. Calculate the propagation factor R
    R_num = T_J * (1 - m_c**2)
    R_den = 1 - (T_J * m_c)**2
    if R_den == 0:
        print("Error: Division by zero in R calculation. Parameters may lead to instability.")
        return
    R = R_num / R_den
    
    # 4. Check for stability condition
    if abs((c - 1) * R) >= 1:
        print("Warning: Stability condition |(c-1)R| < 1 is not met.")
        print(f"|(c-1)R| = {abs((c - 1) * R):.4f}")
    
    # 5. Calculate the constant N
    N = beta * (c * (1 - m0**2)) / (c - 1)

    # 6. Calculate chi
    chi_num = N * (c - 1) * R
    chi_den = 1 - (c - 1) * R
    if chi_den == 0:
        print("Error: Division by zero in chi calculation. The system is at a critical point.")
        return
    chi = chi_num / chi_den

    # 7. Print the results and the final equation with numbers
    print("--- Calculated Parameters ---")
    print(f"T_J = tanh(beta*J) = {T_J:.4f}")
    print(f"Cavity magnetization m_c = {m_c:.4f}")
    print(f"Propagation factor R = {R:.4f}")
    print(f"Constant N = {N:.4f}")
    print("\n--- Final Equation for Susceptibility (chi) ---")
    print("chi = [N * (c-1) * R] / [1 - (c-1) * R]")
    
    # Final output string with formatted numbers for the equation
    final_equation = (
        f"chi = [{N:.4f} * {c-1} * {R:.4f}] / [1 - {c-1} * {R:.4f}] = "
        f"{chi_num:.4f} / {chi_den:.4f} = {chi:.4f}"
    )
    print(final_equation)


if __name__ == '__main__':
    # --- Example Parameters ---
    # Feel free to change these values to explore different scenarios.
    
    # Connectivity of the random graph
    c_param = 3.0
    # Coupling constant
    J_param = 0.5
    # Inverse temperature
    beta_param = 1.0
    # Site magnetization (determined by c, J, beta, but taken as input here)
    m0_param = 0.82
    
    print("--- Input Parameters ---")
    print(f"c = {c_param}, J = {J_param}, beta = {beta_param}, m0 = {m0_param}\n")
    
    calculate_susceptibility(c=c_param, J=J_param, beta=beta_param, m0=m0_param)