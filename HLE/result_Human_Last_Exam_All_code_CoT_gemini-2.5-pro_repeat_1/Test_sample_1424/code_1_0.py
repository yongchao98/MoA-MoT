import cmath

def calculate_bare_greens_function():
    """
    Calculates the retarded bare Green's function G_0 for a given set of parameters.
    """
    # Define parameters for the calculation
    omega = 2.5  # Frequency of the particle (in some energy units)
    epsilon_k = 2.0  # Single-particle energy eigenvalue for state 'k'
    eta = 0.1  # A small positive value for the imaginary part

    # Explain the functional dependence
    print("The retarded bare Green's function G_0(k, omega) depends on the single-particle energy")
    print("eigenvalue epsilon_k according to the formula:")
    print("\n  G_0(k, omega) = 1 / (omega - epsilon_k + i*eta)\n")
    print("This shows an inverse relationship between G_0 and the term (omega - epsilon_k).\n")

    # --- Perform and display a sample calculation ---
    print("--- Example Calculation ---")
    print(f"Given the parameters:")
    print(f"  Frequency (omega) = {omega}")
    print(f"  Energy Eigenvalue (epsilon_k) = {epsilon_k}")
    print(f"  Infinitesimal (eta) = {eta}\n")

    # Calculate the denominator of the Green's function
    denominator = (omega - epsilon_k) + 1j * eta

    # Calculate the Green's function
    G0 = 1 / denominator

    # Print the final equation with all the numbers
    print("The final equation is:")
    print(f"G_0 = 1 / (({omega} - {epsilon_k}) + {eta}j)")
    print(f"G_0 = 1 / ({denominator.real:.2f} + {denominator.imag:.2f}j)")
    print(f"\nThe calculated value of the bare Green's function is:")
    print(f"G_0 = {G0.real:.4f} + ({G0.imag:.4f})j")

if __name__ == "__main__":
    calculate_bare_greens_function()