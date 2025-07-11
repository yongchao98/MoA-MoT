import numpy as np

def calculate_susceptibility(beta, J, c, m0, mcav):
    """
    Calculates the magnetic susceptibility chi for an Ising model on a sparse random graph.

    The final formula derived is:
    chi = N * (c - 1) * K / (1 - (c - 1) * K)
    where:
    N = beta * c * (1 - m0**2) / (c - 1)
    K = tanh(beta * J) * (1 - mcav**2) / (1 - tanh(beta * J)**2 * mcav**2)

    Args:
        beta (float): Inverse temperature (1/kT).
        J (float): Homogeneous coupling constant.
        c (float): Connectivity of the graph.
        m0 (float): Magnetization of the central spin <sigma_0>.
        mcav (float): Cavity magnetization.

    Returns:
        float: The magnetic susceptibility chi, or infinity if the system is unstable.
    """
    if c <= 2:
        raise ValueError("Connectivity c must be greater than 2.")

    # Calculate the propagation factor K
    t = np.tanh(beta * J)
    K_numerator = t * (1 - mcav**2)
    K_denominator = 1 - (t**2) * (mcav**2)
    K = K_numerator / K_denominator

    # The denominator of the susceptibility formula
    chi_denominator = 1 - (c - 1) * K
    
    # Check for the stability condition (paramagnetic phase)
    if chi_denominator <= 0:
        print("Warning: The system is at or beyond the critical point. Susceptibility diverges.")
        return float('inf')

    # Calculate the constant N as defined in the problem
    N = beta * c * (1 - m0**2) / (c - 1)

    # Calculate chi using the final derived formula
    chi_numerator = N * (c-1) * K
    chi = chi_numerator / chi_denominator

    # As requested, output the components of the final equation
    print("--- Calculation Details ---")
    print(f"Final formula: chi = N * (c-1) * K / (1 - (c-1) * K)")
    print("\nComponent values:")
    print(f"N = {N}")
    print(f"K = {K}")
    print(f"c = {c}")
    print(f"Numerator = N * (c - 1) * K = {chi_numerator}")
    print(f"Denominator = 1 - (c - 1) * K = {chi_denominator}")
    print("\n--- Final Result ---")
    print(f"chi = {chi}")
    
    return chi

# Example: Paramagnetic phase (B=0), where m0 and mcav are both 0.
# We use hypothetical values for beta, J, and c for a demonstration.
if __name__ == '__main__':
    # User-definable parameters
    beta_val = 0.5  # Inverse temperature
    J_val = 0.4     # Coupling constant
    c_val = 3.0     # Connectivity
    
    # For this example, we assume we are in the paramagnetic phase where B=0, m0=0, mcav=0
    m0_val = 0.0
    mcav_val = 0.0
    
    calculate_susceptibility(beta=beta_val, J=J_val, c=c_val, m0=m0_val, mcav=mcav_val)
