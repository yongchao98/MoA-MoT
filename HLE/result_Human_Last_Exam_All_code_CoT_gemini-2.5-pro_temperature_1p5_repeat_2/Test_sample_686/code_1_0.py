import numpy as np
from scipy.optimize import fsolve

def solve_ising_susceptibility():
    """
    Calculates the magnetic susceptibility of an Ising model on a sparse random graph.
    """
    # System parameters (user can modify these)
    c = 3.0      # Connectivity
    J = 0.5      # Coupling constant
    B = 0.1      # External field
    beta = 1.0   # Inverse temperature

    t = np.tanh(beta * J)

    # Self-consistency equation for cavity magnetization m_cav
    # f(m_cav) = m_cav - tanh(beta*B + (c-1)*arctanh(t*m_cav)) = 0
    # We must be careful with the domain of arctanh, |t*m_cav| < 1.
    # Since |t|<1 and |m_cav|<1, this is guaranteed.
    def m_cav_equation(m_cav_var):
        m_cav = m_cav_var[0]
        # Use np.arctanh which handles arrays, though we only have a scalar.
        # Ensure the argument is within the valid range to avoid runtime warnings, although mathematically it should be.
        arg = t * m_cav
        if abs(arg) >= 1.0:
            arg = np.sign(arg) * (1.0 - 1e-9)
        return m_cav - np.tanh(beta * B + (c - 1) * np.arctanh(arg))

    # Numerically solve for m_cav
    # Initial guess for m_cav is tanh(beta*B)
    m_cav_initial_guess = [np.tanh(beta*B)]
    m_cav_solution = fsolve(m_cav_equation, m_cav_initial_guess)
    m_cav = m_cav_solution[0]

    # Calculate m0 using the found m_cav
    m0 = (m_cav * (1 + t)) / (1 + t * m_cav**2)

    # Calculate the propagation factor T
    T = t * (1 - m_cav**2) / (1 - (t * m_cav)**2)

    # Calculate susceptibility chi
    denominator = 1 - (c - 1) * T
    if denominator <= 0:
        print("Warning: The parameters lead to a divergent susceptibility (system might be at or below a critical point).")
        chi = float('inf')
    else:
        chi = beta * c * (1 - m0**2) * T / denominator

    # Output the results and the final equation
    print("For the given parameters:")
    print(f"c = {c}, J = {J}, B = {B}, beta = {beta}")
    print("\nIntermediate calculated values:")
    print(f"t = tanh(beta*J) = {t:.4f}")
    print(f"m_cav = {m_cav:.4f}")
    print(f"m0 = {m0:.4f}")
    print(f"T = {T:.4f}")

    print("\nThe final equation for susceptibility is:")
    print(f"chi = beta * c * (1 - m0^2) * T / (1 - (c-1)*T)")
    print("\nSubstituting the numerical values:")
    print(f"chi = {beta:.4f} * {c:.4f} * (1 - {m0:.4f}^2) * {T:.4f} / (1 - ({c-1:.4f})*{T:.4f})")
    print(f"chi = {beta:.4f} * {c:.4f} * ({1-m0**2:.4f}) * {T:.4f} / ({denominator:.4f})")
    
    numerator_val = beta * c * (1 - m0**2) * T
    
    print(f"chi = {numerator_val:.4f} / {denominator:.4f}")
    
    print(f"\nFinal result:")
    print(f"chi = {chi:.4f}")
    
    # Return the final result in the requested format
    print(f"\n<<<{chi:.4f}>>>")

if __name__ == '__main__':
    solve_ising_susceptibility()
