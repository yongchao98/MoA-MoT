import cmath

def calculate_bare_greens_function(i_omega_n, epsilon_k):
    """
    Calculates the bare Green's function in Matsubara frequency space.

    Args:
        i_omega_n (complex): The fermionic Matsubara frequency (e.g., 3.0j).
        epsilon_k (float): The single-particle energy eigenvalue.

    Returns:
        complex: The value of the bare Green's function G_0.
    """
    # The core formula for the bare Green's function
    G_0 = 1 / (i_omega_n - epsilon_k)
    return G_0

def main():
    """
    Main function to demonstrate the calculation and print the results.
    """
    # --- Parameters ---
    # Example single-particle energy eigenvalue (relative to the chemical potential)
    # Units are arbitrary (e.g., eV).
    epsilon_k = 2.0

    # Example fermionic Matsubara frequency. In Python, the imaginary unit is 'j'.
    # This would correspond to a specific temperature and integer n in omega_n = (2n+1)pi*T
    i_omega_n = 3.0j

    # --- Calculation ---
    G0_value = calculate_bare_greens_function(i_omega_n, epsilon_k)

    # --- Output ---
    print("In the Matsubara formalism, the bare Green's function G_0 has a direct dependence on the single-particle energy epsilon_k.")
    print("The functional form is: G_0(k, i*omega_n) = 1 / (i*omega_n - epsilon_k)")
    print("\nLet's calculate an example:")

    # Print the equation with the specific numbers used
    print(f"\nFor epsilon_k = {epsilon_k} and i*omega_n = {i_omega_n}:")
    print(f"G_0 = 1 / ({i_omega_n} - {epsilon_k})")

    # Print the final calculated value
    print(f"\nG_0 = {G0_value}")

if __name__ == "__main__":
    main()
