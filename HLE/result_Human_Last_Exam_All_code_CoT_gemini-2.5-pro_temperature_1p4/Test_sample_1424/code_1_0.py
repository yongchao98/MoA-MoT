import cmath

def demonstrate_green_function(omega, epsilon_k, eta=0.1):
    """
    Calculates the bare Green's function G_0 for given parameters
    and prints the functional relationship.

    The bare Green's function (retarded) is defined as:
    G_0(k, ω) = 1 / (ω - ε_k + iη)

    Args:
        omega (float): The frequency.
        epsilon_k (float): The single-particle energy eigenvalue.
        eta (float): A small positive infinitesimal for calculation.
    """
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy ε_k is given by:")
    print("G_0(k, ω) = 1 / (ω - ε_k + iη)\n")

    print("Let's demonstrate with a numerical example:")
    print(f"  Frequency (ω) = {omega}")
    print(f"  Energy Eigenvalue (ε_k) = {epsilon_k}")
    print(f"  Infinitesimal (η) = {eta}\n")
    
    # Calculate the denominator and the final G_0 value
    denominator = (omega - epsilon_k) + 1j * eta
    G0_value = 1 / denominator

    print("Substituting these values into the equation:")
    # Here we output each number in the final equation, as requested
    print(f"  G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")
    print(f"  G_0 = 1 / {denominator}")
    print(f"  G_0 = {G0_value}")
    print(f"  G_0 ≈ {G0_value.real:.4f} + {G0_value.imag:.4f}i")


# --- Main Execution ---
# You can change these values to see how G_0 changes
example_omega = 2.5
example_epsilon_k = 2.0

demonstrate_green_function(example_omega, example_epsilon_k)