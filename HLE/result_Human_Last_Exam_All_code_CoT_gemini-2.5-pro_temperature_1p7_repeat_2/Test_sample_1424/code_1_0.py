import cmath

def calculate_bare_greens_function():
    """
    Calculates and displays the bare Green's function G_0 for a given set of parameters.
    """
    # --- Parameters ---
    # In quantum mechanics, energy and frequency are related by E = ħω.
    # For simplicity in many theoretical calculations, it's common to use units where ħ=1.
    hbar = 1.0

    # Let's choose an example frequency and a single-particle energy eigenvalue.
    # omega represents the frequency (related to the energy of the propagating particle).
    omega = 2.5
    # epsilon_k is the energy of a single-particle state with momentum 'k'.
    epsilon_k = 2.0
    # delta is a small positive infinitesimal. It is crucial for defining the contour
    # of integration in the complex plane and ensures causality.
    delta = 0.1

    # --- The Formula ---
    # The Feynman Green's function G_0 has a specific dependence on epsilon_k.
    # G_0(k, ω) = 1 / (ħω - ε_k + iδ)
    # The core dependence is an inverse relationship with the energy difference (ħω - ε_k).

    print("The bare Green's function G_0(k, ω) is given by the formula:")
    print("  G_0 = 1 / (ħω - ε_k + iδ)\n")
    print("This shows G_0 is inversely proportional to the term (ħω - ε_k).\n")
    print("--- Example Calculation ---")
    print(f"Using the following values:")
    print(f"  ħ       = {hbar}")
    print(f"  ω       = {omega}")
    print(f"  ε_k     = {epsilon_k}")
    print(f"  δ       = {delta}")
    print(f"  i (j)   = {1j}\n")

    # --- Calculation Steps ---
    print("Substituting the values into the equation:")
    # First, calculate the term in the denominator
    energy_term = hbar * omega
    denominator = (energy_term - epsilon_k + 1j * delta)
    
    # Calculate the final Green's function
    G0 = 1 / denominator

    # Print the step-by-step equation with numbers
    print(f"G_0 = 1 / (({hbar} * {omega}) - {epsilon_k} + {delta}j)")
    print(f"G_0 = 1 / ({energy_term} - {epsilon_k} + {delta}j)")
    print(f"G_0 = 1 / ({denominator.real} + {denominator.imag}j)")

    # Print the final result in a standard complex number format
    result_str = f"{G0.real:.4f} - {-G0.imag:.4f}j" if G0.imag < 0 else f"{G0.real:.4f} + {G0.imag:.4f}j"
    print(f"\nThe calculated value for G_0 is: {result_str}")


if __name__ == "__main__":
    calculate_bare_greens_function()