import cmath

def calculate_bare_greens_function(omega, epsilon_k, delta):
    """
    Calculates the bare Green's function G_0 in frequency-momentum space.

    Args:
        omega (float): The frequency/energy of the particle.
        epsilon_k (float): The single-particle energy eigenvalue.
        delta (float): A small positive infinitesimal for contour integration.

    Returns:
        complex: The value of the bare Green's function G_0.
    """
    # The Green's function is G_0 = 1 / (omega - epsilon_k + i*delta)
    denominator = (omega - epsilon_k) + 1j * delta
    g0 = 1 / denominator
    return g0

def main():
    """
    Main function to demonstrate the calculation of G_0.
    """
    # --- Example values ---
    # Frequency/energy of the propagating particle
    omega = 2.5
    # Single-particle energy of a specific state 'k'
    epsilon_k = 2.4
    # Small infinitesimal value
    delta = 0.1

    # Calculate G_0
    g0_value = calculate_bare_greens_function(omega, epsilon_k, delta)

    # Print the equation with the numbers used
    print("The bare Green's function G_0 is calculated using the formula:")
    print("G_0 = 1 / (omega - epsilon_k + i*delta)")
    print("\nWith the following values:")
    print(f"omega     (ω) = {omega}")
    print(f"epsilon_k (ε_k) = {epsilon_k}")
    print(f"delta     (δ) = {delta}")
    
    # Final equation with numbers
    print("\nFinal Equation:")
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{delta})")

    # Print the result
    print(f"\nThe calculated value of G_0 is: {g0_value}")
    
    # You can also print the real and imaginary parts separately
    print(f"Real part: {g0_value.real:.4f}, Imaginary part: {g0_value.imag:.4f}")

if __name__ == "__main__":
    main()
