import numpy as np
from scipy.optimize import fsolve

def calculate_susceptibility(c, J, beta, B):
    """
    Calculates the magnetic susceptibility of an Ising model on a sparse random graph.

    Args:
        c (float): Connectivity of the graph.
        J (float): Coupling constant.
        beta (float): Inverse temperature (1/kT).
        B (float): External magnetic field.
    """
    # 1. Solve for the self-consistent message 'u'
    # The equation is u = (1/beta) * arctanh[tanh(beta*J) * tanh(beta*(B + (c-1)*u))]
    # We find the root of f(u) = 0.
    def message_equation(u):
        # Use np.clip to avoid tanh(>1) issues in arctanh due to numerical precision
        tanh_arg = np.tanh(beta * (B + (c - 1) * u))
        inner_term = np.clip(np.tanh(beta * J) * tanh_arg, -0.99999999, 0.99999999)
        return u - (1 / beta) * np.arctanh(inner_term)

    # Initial guess for u. 0 is a good guess for small B.
    u_initial_guess = 0.0
    u_solution = fsolve(message_equation, u_initial_guess)[0]

    # 2. Calculate magnetizations and the propagation factor T
    m_cav = np.tanh(beta * (B + (c - 1) * u_solution))
    m0 = np.tanh(beta * (B + c * u_solution))
    T = beta * (1 - m_cav**2) * np.tanh(beta * J)

    # 3. Calculate susceptibility chi
    numerator = c * beta * (1 - m0**2) * T
    denominator = 1 - (c - 1) * T

    if abs(denominator) < 1e-9:
        print("Warning: Denominator is close to zero. System may be at or near a critical point.")
        chi = float('inf')
    else:
        chi = numerator / denominator

    # 4. Print the components of the final equation and the result
    print("--- Calculation Parameters ---")
    print(f"Connectivity c = {c}")
    print(f"Coupling J = {J}")
    print(f"Inverse Temperature beta = {beta}")
    print(f"External Field B = {B}")
    print("\n--- Intermediate Values ---")
    print(f"Self-consistent message u = {u_solution}")
    print(f"Cavity magnetization m_cav = {m_cav}")
    print(f"Site magnetization m_0 = {m0}")
    print(f"Propagation factor T = {T}")
    print("\n--- Final Equation Components ---")
    print(f"χ = (c * β * (1 - m_0²) * T) / (1 - (c - 1) * T)")
    print(f"Numerator = {c} * {beta:.4f} * (1 - {m0:.4f}²) * {T:.4f} = {numerator:.4f}")
    print(f"Denominator = 1 - ({c-1}) * {T:.4f} = {denominator:.4f}")
    print("\n--- Result ---")
    print(f"Magnetic Susceptibility χ = {chi}")


if __name__ == '__main__':
    # --- Example Parameters ---
    # Connectivity (must be > 2)
    c_param = 3.0
    # Coupling constant
    J_param = 1.0
    # Inverse temperature (beta = 1/kT).
    # For c=3, J=1, the critical point is at beta_crit = arctanh(1/(c-1))/J = arctanh(0.5) approx 0.549
    beta_param = 0.5
    # External field (set to 0 for zero-field susceptibility)
    B_param = 0.0

    calculate_susceptibility(c=c_param, J=J_param, beta=beta_param, B=B_param)