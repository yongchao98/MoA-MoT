import cmath

def calculate_bare_greens_function():
    """
    Calculates the bare Green's function for a given set of parameters
    and prints the calculation step-by-step.
    """
    # Define the parameters for the calculation
    # omega: the energy of the probe
    omega = 1.5  # in some arbitrary units of energy, e.g., eV
    
    # epsilon_k: the energy eigenvalue of the single-particle state 'k'
    epsilon_k = 1.2 # in the same units as omega
    
    # eta: a small positive number for the infinitesimal imaginary part
    eta = 0.1

    # The formula is G_0 = 1 / (omega - epsilon_k + i*eta)
    print(f"The functional form of the bare Green's function is G_0 = 1 / (omega - epsilon_k + i*eta)")
    print("-" * 30)
    print("Calculating for the following values:")
    print(f"omega     = {omega}")
    print(f"epsilon_k = {epsilon_k}")
    print(f"eta       = {eta}")
    print("-" * 30)

    # Calculate the denominator
    denominator = (omega - epsilon_k) + 1j * eta
    
    # Calculate the Green's function
    G0 = 1 / denominator

    # Print the final equation with all the numbers
    print("The final equation is:")
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")
    print(f"G_0 = 1 / {denominator}")
    print(f"G_0 = {G0}")

if __name__ == "__main__":
    calculate_bare_greens_function()