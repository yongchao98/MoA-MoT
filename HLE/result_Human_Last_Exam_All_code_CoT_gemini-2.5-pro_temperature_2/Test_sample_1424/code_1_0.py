import cmath

def calculate_bare_greens_function():
    """
    Calculates the bare Green's function G_0 for given parameters
    and prints the calculation step-by-step.
    """
    # --- Parameters ---
    # omega: The frequency (energy) of the probe
    # epsilon_k: The single-particle energy eigenvalue
    # mu: The chemical potential
    # eta: A small positive constant for causality
    omega = 2.5
    epsilon_k = 2.0
    mu = 0.5
    eta = 0.1

    # --- Calculation ---
    # The formula for the retarded bare Green's function is:
    # G_0(k, omega) = 1 / (omega - (epsilon_k - mu) + i*eta)

    # Calculate the energy relative to the chemical potential
    energy_relative = epsilon_k - mu
    
    # Calculate the denominator of the Green's function
    denominator = (omega - energy_relative) + 1j * eta

    # Calculate the Green's function
    G0 = 1 / denominator

    # --- Output ---
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy epsilon_k is:")
    print("G_0(k, omega) = 1 / (omega - (epsilon_k - mu) + i*eta)\n")
    
    print(f"Given the parameters:")
    print(f"  omega     = {omega}")
    print(f"  epsilon_k = {epsilon_k}")
    print(f"  mu        = {mu}")
    print(f"  eta       = {eta}\n")

    print("The calculation is:")
    print(f"G_0 = 1 / ( {omega} - ( {epsilon_k} - {mu} ) + i * {eta} )")
    print(f"G_0 = 1 / ( {omega} - {energy_relative} + i * {eta} )")
    print(f"G_0 = 1 / ( {omega - energy_relative} + i * {eta} )")
    print(f"G_0 = {G0}")

calculate_bare_greens_function()