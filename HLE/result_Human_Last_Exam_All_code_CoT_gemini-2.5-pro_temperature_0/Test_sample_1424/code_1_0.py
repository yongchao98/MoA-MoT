import numpy as np

def calculate_bare_greens_function():
    """
    Calculates and displays the bare Green's function G_0 for a given set of parameters.
    """
    # --- Parameters ---
    # omega: The frequency (energy) of the particle being propagated.
    # epsilon_k: The energy eigenvalue of a specific single-particle state 'k'.
    # mu: The chemical potential, which defines the Fermi level (energy separating occupied and unoccupied states).
    # eta: A small positive infinitesimal to enforce causality.
    omega = 2.5
    epsilon_k = 2.0
    mu = 0.0
    eta = 0.1

    # --- Explanation ---
    print("The bare (non-interacting) time-ordered Green's function G_0(omega, k) depends on the")
    print("single-particle energy eigenvalue epsilon_k as follows:")
    print("\n  G_0(omega, k) = 1 / (omega - epsilon_k + i * eta * sgn(epsilon_k - mu))\n")
    print("This shows G_0 is inversely proportional to the difference between the frequency 'omega'")
    print("and the single-particle energy 'epsilon_k'.\n")
    print("-" * 60)

    # --- Calculation for a specific case ---
    print("Let's calculate G_0 for a specific set of values:")
    print(f"  Frequency (omega)              = {omega}")
    print(f"  Single-particle energy (eps_k) = {epsilon_k}")
    print(f"  Chemical potential (mu)        = {mu}")
    print(f"  Infinitesimal (eta)            = {eta}\n")

    # Calculate the sign term, which determines if the state is a particle or hole
    sign_term = np.sign(epsilon_k - mu)

    # --- Output the final equation with numbers ---
    # Determine the sign for the imaginary part in the string representation
    imag_sign = "+" if sign_term >= 0 else "-"
    # The value of the imaginary part in the denominator
    imag_value = eta * abs(sign_term)

    print("The final equation with these numbers is:")
    # Using 1j for the imaginary unit in the calculation
    denominator_calc = (omega - epsilon_k) + 1j * eta * sign_term
    final_result = 1 / denominator_calc

    # We explicitly show each number in the printed equation
    print(f"  G_0 = 1 / ( {omega} - {epsilon_k} {imag_sign} i * {imag_value} )")
    print(f"  G_0 = 1 / ( {omega - epsilon_k} {imag_sign} i * {imag_value} )")
    print("\n" + "-" * 60)
    print("The final calculated value is:")
    print(f"  G_0 = {final_result}")

if __name__ == "__main__":
    calculate_bare_greens_function()