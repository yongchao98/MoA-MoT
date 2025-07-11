import numpy as np

def calculate_bare_greens_function():
    """
    Calculates and demonstrates the bare Green's function G_0.
    """
    # --- Define example parameters ---
    omega = 2.5       # Frequency
    epsilon_k = 2.0   # Single-particle energy eigenvalue for state k
    mu = 0.5          # Chemical potential
    delta = 1e-9      # A small positive infinitesimal for causality

    # --- Print explanations ---
    print("The functional dependence of the bare Green's function G_0 on the single-particle")
    print("energy eigenvalue (epsilon_k) is given by the formula:")
    print("\n  G_0(k, omega) = 1 / (omega - epsilon_k + i*delta*sgn(epsilon_k - mu))\n")
    print("This shows G_0 has a simple pole where the real part is determined by epsilon_k.\n")

    # --- Calculation for the example case ---
    print("--- Calculating for an example case ---")
    print(f"Frequency (omega)                  = {omega}")
    print(f"Single-particle energy (epsilon_k) = {epsilon_k}")
    print(f"Chemical Potential (mu)            = {mu}")
    print(f"Infinitesimal (delta)              = {delta}\n")

    # Calculate the sign term
    # This determines if we are dealing with a particle (positive) or hole (negative)
    sign_term = np.sign(epsilon_k - mu)

    # --- Display the final equation with all numbers plugged in ---
    print("Substituting these values into the formula step-by-step:")

    # Format the imaginary part's sign for clear output
    imaginary_part_sign_str = "+" if sign_term >= 0 else "-"
    
    # Calculate the final denominator and the Green's function
    denominator = (omega - epsilon_k) + 1j * delta * sign_term
    G0 = 1 / denominator

    # Print the equation with all the variables from the problem statement
    print(f"G_0 = 1 / ( omega - epsilon_k + i * delta * sgn(epsilon_k - mu) )")
    print(f"G_0 = 1 / ( {omega} - {epsilon_k} + i * {delta} * sgn({epsilon_k} - {mu}) )")
    print(f"G_0 = 1 / ( {omega} - {epsilon_k} + i * {delta} * ({int(sign_term)}) )")
    print(f"G_0 = 1 / ( {omega - epsilon_k} {imaginary_part_sign_str} i * {delta} )")
    print(f"G_0 = 1 / ( {denominator} )")
    print(f"\nFinal Result: G_0 = {G0}")

# Run the calculation and print the results
calculate_bare_greens_function()
