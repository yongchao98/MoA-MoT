import math

def evaluate_conductivity_correction():
    """
    Calculates the quantum correction to conductivity for an electron in a 3D
    disordered medium (weak localization).
    """
    # Physical constants
    e = 1.602176634e-19  # Elementary charge in Coulombs
    hbar = 1.054571817e-34 # Reduced Planck constant in J·s
    pi = math.pi

    # Typical parameters for a doped semiconductor at low temperature
    # l_e: elastic mean free path in meters
    l_e = 50e-9
    # L_phi: phase coherence length in meters
    # This must be greater than l_e for the theory to apply.
    L_phi = 500e-9

    print("Evaluating the quantum correction to conductivity (Weak Localization in 3D):")
    print("Formula: Δσ = - (e² / (2 * π² * ħ)) * (1/l_e - 1/L_φ)\n")

    # Print the formula with numerical values
    print("Plugging in the values:")
    print(f"Δσ = - (({e:.4e} C)² / (2 * π² * {hbar:.4e} J·s)) * (1 / {l_e:.2e} m - 1 / {L_phi:.2e} m)\n")

    # Calculate each part of the equation
    prefactor_val = e**2 / (2 * pi**2 * hbar)
    length_term_val = (1/l_e - 1/L_phi)
    
    # Print the intermediate calculated values
    print("Calculating the prefactor and the length-dependent term:")
    print(f"Prefactor (e² / (2 * π² * ħ)) = {prefactor_val:.4e} S (Siemens)")
    print(f"Length term (1/l_e - 1/L_φ) = ({1/l_e:.2e} - {1/L_phi:.2e}) m⁻¹ = {length_term_val:.2e} m⁻¹\n")

    # Calculate the final result
    delta_sigma = -prefactor_val * length_term_val

    # Print the final calculation step and result
    print("Final Calculation:")
    print(f"Δσ = - ({prefactor_val:.4e} S) * ({length_term_val:.2e} m⁻¹)")
    print(f"Δσ = {delta_sigma:.2f} S/m")
    
    return delta_sigma

if __name__ == '__main__':
    final_answer = evaluate_conductivity_correction()
    # The final answer is returned in the specified format below.
    # <<<final_answer>>>
    
<<<%.2f>>> % final_answer