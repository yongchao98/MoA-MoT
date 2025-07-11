import math

def evaluate_quantum_correction_formula():
    """
    This function prints the formula for the quantum correction to conductivity
    in a bulk (3D) semiconductor, also known as the weak localization correction.
    """

    # Define the symbols for the formula using unicode for better representation
    e = "e"
    pi = "π"
    h_bar = "ħ"
    l = "l"
    L_phi = "L_φ"
    delta_sigma = "δσ"

    # Construct the formula string. The formula contains the number '2' as requested.
    formula_prefactor = f"-({e}² / (2 * {pi}² * {h_bar}))"
    formula_lengths = f"(1/{l} - 1/{L_phi})"
    full_formula = f"{delta_sigma} = {formula_prefactor} * {formula_lengths}"

    # Print the explanation and the formula
    print("The quantum correction to conductivity (weak localization) for a bulk (3D) system is evaluated by the following formula:")
    print("\n" + full_formula + "\n")

    # Print the definitions of the symbols used in the formula
    print("Where:")
    print(f"  {delta_sigma} : The quantum correction to conductivity.")
    print(f"  {e} : The elementary charge of an electron.")
    print(f"  {pi} : The mathematical constant Pi ({math.pi:.4f}...).")
    print(f"  {h_bar} : The reduced Planck's constant (h/2π).")
    print(f"  {l} : The elastic mean free path, which acts as the short-distance cutoff for interference.")
    print(f"  {L_phi} : The phase coherence length, the long-distance cutoff beyond which interference effects are suppressed.")

# Execute the function to display the result
evaluate_quantum_correction_formula()
