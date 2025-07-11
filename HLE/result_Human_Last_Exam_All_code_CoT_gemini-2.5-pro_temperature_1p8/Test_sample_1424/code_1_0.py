import cmath

def display_bare_greens_function(omega, epsilon_k, delta=0.1):
    """
    Calculates and displays the bare Green's function G_0(k, omega).

    In the Feynman path integral formalism, the bare Green's function G_0 has a
    simple inverse dependence on the single-particle energy eigenvalue epsilon_k.
    The formula in frequency-momentum space is:

        G_0(k, omega) = 1 / (omega - epsilon_k + i*delta)

    where delta is a small positive infinitesimal that enforces causality.

    This function demonstrates this relationship with example values.

    Args:
        omega (float): The frequency.
        epsilon_k (float): The single-particle energy eigenvalue.
        delta (float, optional): A small positive value for the imaginary part.
                                 Defaults to 0.1 for visualization.
    """
    # The term in the denominator
    # For calculation, we use the cmath library for complex numbers
    denominator = (omega - epsilon_k) + 1j * delta
    
    # Calculate the Green's function
    G_0 = 1 / denominator

    # Print the explanation and the final equation with all numbers
    print("The bare Green's function G_0 is a function of frequency (omega) and the single-particle energy eigenvalue (epsilon_k).")
    print("The functional dependence is an inverse relationship:")
    print("G_0(k, omega) = 1 / (omega - epsilon_k + i*delta)\n")
    print("For the given values:")
    print(f"  omega     = {omega}")
    print(f"  epsilon_k = {epsilon_k}")
    print(f"  delta     = {delta}\n")
    print("The equation is:")
    # We use .real and .imag for clean printing of the complex denominator
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + {denominator.imag:.1f}i)")
    print(f"G_0 = 1 / ({denominator.real} + {denominator.imag:.1f}i)")
    print(f"G_0 = {G_0.real:.4f} + {G_0.imag:.4f}i")


if __name__ == '__main__':
    # --- Example Values ---
    # Frequency of the particle
    example_omega = 5.0
    # Energy of the single-particle state 'k'
    example_epsilon_k = 3.5

    # Calculate and display the result
    display_bare_greens_function(example_omega, example_epsilon_k)
