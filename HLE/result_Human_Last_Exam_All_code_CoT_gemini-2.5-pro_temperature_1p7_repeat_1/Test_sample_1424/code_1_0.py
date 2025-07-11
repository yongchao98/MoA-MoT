import cmath

def calculate_bare_greens_function(omega, epsilon_k, eta=0.1):
    """
    Calculates the bare Green's function G_0 in frequency-momentum space.

    Args:
        omega (float): The frequency variable.
        epsilon_k (float): The single-particle energy eigenvalue.
        eta (float): A small positive infinitesimal for the Feynman prescription.
    """
    # The functional dependence is G_0 = 1 / (ω - ε_k + iη)
    denominator = (omega - epsilon_k) + 1j * eta
    g0 = 1 / denominator

    # Print the explanation and the step-by-step calculation
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy ε_k is given by:")
    print("G_0(k, ω) = 1 / (ω - ε_k + iη)")
    print("\nwhere:")
    print("  ω: frequency")
    print("  ε_k: single-particle energy eigenvalue for momentum k")
    print("  iη: a small imaginary part to ensure causality (Feynman prescription)")
    print("-" * 30)

    # Outputting each number in the final equation as requested
    print("Let's calculate for a sample set of values:")
    print(f"ω = {omega}")
    print(f"ε_k = {epsilon_k}")
    print(f"η = {eta}")
    print("\nPlugging the numbers into the equation:")
    # The final equation with each number explicitly shown
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")
    print(f"G_0 = 1 / {denominator}")
    print(f"G_0 = {g0}")

# --- Example Usage ---
# Define sample values for the variables
# Let's choose a frequency 'omega' close to the particle energy 'epsilon_k'
# to see a significant response.
omega_val = 2.5
epsilon_k_val = 2.4

# Calculate and display the result
calculate_bare_greens_function(omega_val, epsilon_k_val)
