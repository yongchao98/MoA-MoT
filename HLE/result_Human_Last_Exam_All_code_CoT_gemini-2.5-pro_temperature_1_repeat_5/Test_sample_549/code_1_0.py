import math

def evaluate_quantum_correction():
    """
    This function evaluates and prints the formula for the quantum correction
    to conductivity in a bulk (3D) semiconductor based on weak localization theory.
    """
    print("The quantum correction to conductivity (Δσ) for an electron in a bulk semiconductor is a result of the weak localization effect.")
    print("It is determined by the interplay of elastic scattering and phase-breaking events.")
    print("\nThe established formula for this correction in a 3D system is presented below.")
    print("-" * 50)

    # The formula is: Δσ = - (e^2 / (2 * π^2 * ħ)) * (1/l_e - 1/L_φ)
    # We will print this formula and explain each term.

    # Components of the equation
    term1 = "e²"
    term2 = f"2 * π² * ħ"
    term3 = "1 / l_e"
    term4 = "1 / L_φ"
    
    print("\nFinal Equation:")
    print(f"Δσ = - ( {term1} / ({term2}) ) * ( ({term3}) - ({term4}) )")

    print("\nWhere the terms represent:")
    print("  Δσ: The quantum correction to conductivity.")
    print(f"  {term1}: The square of the elementary charge 'e'.")
    print(f"  2: The number two.")
    print(f"  π²: The square of the mathematical constant Pi ({round(math.pi, 4)}...).")
    print(f"  ħ: The reduced Planck constant (h/2π).")
    print(f"  {term3}: The inverse of the elastic mean free path 'l_e'.")
    print(f"  {term4}: The inverse of the phase coherence length 'L_φ'.")

    print("-" * 50)
    print("\nNote: Since the phase coherence length L_φ must be greater than the elastic mean free path l_e, the term (1/l_e - 1/L_φ) is positive.")
    print("The overall negative sign indicates that weak localization leads to a reduction in conductivity.")

# Execute the function to print the evaluation
evaluate_quantum_correction()
