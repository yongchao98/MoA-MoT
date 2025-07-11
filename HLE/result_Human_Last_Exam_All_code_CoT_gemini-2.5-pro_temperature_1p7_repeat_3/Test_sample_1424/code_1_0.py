import cmath

def get_bare_greens_function(omega, epsilon_k):
    """
    Calculates the bare Green's function G_0 for a given frequency and energy.

    The formula is G_0(k, ω) = 1 / (ω - ϵ_k + iη), where η is an
    infinitesimal positive constant.
    """
    # Define an infinitesimally small positive constant for eta
    eta = 1e-9

    # Calculate the bare Green's function
    # The term 1j is the Python representation for the imaginary unit i
    G0 = 1 / (omega - epsilon_k + 1j * eta)

    # Print the explanation and the final equation with all numbers
    print("The bare Green's function G_0 is calculated using the formula:")
    print("G_0(k, ω) = 1 / (ω - ϵ_k + iη)\n")
    print("For the given values:")
    print(f"ω = {omega}")
    print(f"ϵ_k = {epsilon_k}")
    print(f"η = {eta}\n")
    print("The final equation is:")
    # Using .real and .imag for a clearer representation of the complex number result
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + {eta}i) = {G0.real:.4f} + {G0.imag:.4f}i")

# --- Example ---
# Set example values for the frequency (ω) and the single-particle energy (ϵ_k)
frequency_omega = 3.5
energy_epsilon_k = 3.2

# Run the calculation and print the result
get_bare_greens_function(frequency_omega, energy_epsilon_k)