import math

def evaluate_quantum_correction():
    """
    This function evaluates and prints the formula for the quantum correction
    to conductivity in a bulk (3D) semiconductor based on weak localization theory.
    """

    # Define the symbols used in the formula for clear printing
    delta_sigma = "Δσ"
    e_sq = "e^2"
    hbar = "ħ"
    l_e = "l_e"
    L_phi = "L_φ"
    
    # The formula is Δσ = - (e^2 / (2 * π^2 * ħ)) * (1/l_e - 1/L_φ)
    # The prompt requires outputting each number in the final equation.
    
    print("The formula for the quantum correction to conductivity (Δσ) in a 3D semiconductor is:")
    print(f"{delta_sigma} = - ( {e_sq} / (2 * π^2 * {hbar}) ) * ( 1/{l_e} - 1/{L_phi} )")
    
    print("\nWhere the terms are:")
    print(f"  {delta_sigma}: The quantum correction to conductivity.")
    print(f"  e: The elementary charge.")
    print(f"  {hbar}: The reduced Planck constant.")
    print(f"  π: The mathematical constant Pi.")
    print(f"  {l_e}: The electron's elastic mean free path.")
    print(f"  {L_phi}: The electron's phase coherence length.")

    print("\nEvaluating the numerical constants in the equation's prefactor:")
    
    # Explicitly state the numbers as requested
    num_2 = 2
    num_pi = math.pi
    num_pi_sq = num_pi**2
    denominator_numerical = num_2 * num_pi_sq
    
    print(f"  The number '2' is a constant from the theoretical derivation.")
    print(f"  The number 'π' is approximately {num_pi:.5f}.")
    print(f"  Therefore, 'π^2' is approximately {num_pi_sq:.5f}.")
    print(f"  The complete numerical part of the denominator is 2 * π^2 ≈ {denominator_numerical:.5f}.")
    
    prefactor_numerical = -1 / denominator_numerical
    
    print("\nThe final equation, with the numerical prefactor calculated, is:")
    print(f"{delta_sigma} = {prefactor_numerical:.5f} * ( {e_sq} / {hbar} ) * ( 1/{l_e} - 1/{L_phi} )")

if __name__ == '__main__':
    evaluate_quantum_correction()