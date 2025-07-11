import cmath

def calculate_bare_greens_function(omega, epsilon_k, delta=0.1):
    """
    Calculates the bare Green's function G_0(k, omega).

    Args:
        omega (float): The frequency.
        epsilon_k (float): The single-particle energy eigenvalue.
        delta (float): A small positive constant for causality.

    Returns:
        complex: The value of the Green's function.
    """
    return 1 / (omega - epsilon_k + 1j * delta)

def main():
    """
    Main function to demonstrate the calculation.
    """
    # Example values
    omega = 2.5  # Example frequency in some energy units
    epsilon_k = 2.0  # Example single-particle energy in some energy units
    delta = 0.1  # Small positive constant

    # The functional dependence
    print("In the frequency-momentum domain, the bare Green's function G_0 has the following dependence on the single-particle energy eigenvalue ε_k:")
    print("G_0(ω, ε_k) = 1 / (ω - ε_k + iδ)\n")

    # Perform the calculation
    result = calculate_bare_greens_function(omega, epsilon_k, delta)

    # Print the step-by-step calculation
    print(f"Let's calculate G_0 for the following values:")
    print(f"  Frequency ω = {omega}")
    print(f"  Energy Eigenvalue ε_k = {epsilon_k}")
    print(f"  Infinitesimal term δ = {delta}\n")

    print("Substituting the values into the equation:")
    # Using f-string formatting to display the equation with numbers
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i * {delta})")
    
    # Show the intermediate step
    real_part = omega - epsilon_k
    imag_part_str = f"+ i * {delta}"
    print(f"G_0 = 1 / ({real_part} {imag_part_str})")

    # Show the final result
    print(f"\nThe final value of the bare Green's function is: {result}")

if __name__ == "__main__":
    main()