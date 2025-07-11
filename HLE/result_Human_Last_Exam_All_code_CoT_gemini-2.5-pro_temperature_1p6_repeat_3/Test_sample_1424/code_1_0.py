import cmath

def calculate_bare_greens_function(omega, epsilon_k, eta=0.1):
    """
    Calculates the bare Green's function G_0(k, omega).

    Args:
        omega (float): The frequency (energy) of the probe.
        epsilon_k (float): The single-particle energy eigenvalue.
        eta (float): A small positive infinitesimal for causality.

    Returns:
        complex: The value of the Green's function.
    """
    # The formula is G_0 = 1 / (omega - epsilon_k + i*eta)
    denominator = (omega - epsilon_k) + 1j * eta
    G0 = 1 / denominator
    return G0

def main():
    """
    Main function to demonstrate the calculation for various epsilon_k.
    """
    # Define a fixed frequency and a small eta
    omega = 2.0
    eta = 0.1

    # A list of single-particle energy eigenvalues to test
    epsilon_k_values = [1.0, 1.8, 2.0, 2.2, 3.0]

    print(f"Calculating the bare Green's function G_0 for a fixed frequency omega = {omega}\n")

    # Iterate through each epsilon_k and show the calculation
    for epsilon_k in epsilon_k_values:
        G0 = calculate_bare_greens_function(omega, epsilon_k, eta)
        
        # We print the numbers that form the equation to show the functional dependence
        print(f"For epsilon_k = {epsilon_k:.2f}:")
        print(f"  G_0 = 1 / ({omega:.2f} - {epsilon_k:.2f} + i*{eta:.2f})")
        print(f"  G_0 = {G0.real:.4f} + {G0.imag:.4f}i")
        # Notice the magnitude of G_0 is largest when omega is closest to epsilon_k
        print(f"  |G_0| = {abs(G0):.4f}\n")

if __name__ == "__main__":
    main()