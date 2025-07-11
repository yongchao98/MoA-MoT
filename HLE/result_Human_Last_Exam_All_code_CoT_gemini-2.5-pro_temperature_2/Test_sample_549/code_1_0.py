import math

def evaluate_quantum_correction_conductivity():
    """
    This function evaluates and prints the formula for the quantum correction
    to conductivity (weak localization) in a 3D system.
    """

    # --- Define the symbols used in the equation ---
    e = "e"         # Elementary charge
    hbar = "ħ"      # Reduced Planck Constant
    pi = "π"        # The constant Pi
    l_e = "l_e"     # Elastic mean free path
    L_phi = "L_p"   # Phase coherence length (using 'p' for 'phi')

    # --- Construct the full equation string ---
    # The formula is: Δσ = (e² / (π² * ħ)) * [ (1/L_p) * arctan(L_p / l_e) - (1/l_e) ]
    
    # Pre-factor
    factor_numerator = f"{e}²"
    factor_denominator = f"({pi}² * {hbar})"
    factor = f"({factor_numerator} / {factor_denominator})"

    # Term inside the brackets
    term_a = f"(1/{L_phi}) * arctan({L_phi} / {l_e})"
    term_b = f"1/{l_e}"
    bracket_expression = f"[{term_a} - {term_b}]"

    final_equation = f"Δσ = {factor} * {bracket_expression}"

    print("The quantum correction to conductivity, Δσ, for an electron in a bulk (3D) semiconductor is given by the following equation:")
    print("\n" + final_equation + "\n")

    print("Where each component in the equation represents:")
    print(f"Δσ: The quantum correction to the conductivity")
    print(f"{e}: The elementary charge ({1.602e-19:.3e} Coulombs)")
    print(f"{pi}: The mathematical constant Pi ({math.pi:.4f})")
    print(f"{hbar}: The reduced Planck constant ({1.054e-34:.3e} J·s)")
    print(f"{l_e}: The electron's elastic mean free path (the average distance an electron travels between scattering events)")
    print(f"{L_phi}: The electron's phase coherence length (the distance over which the electron's wave function maintains its phase before decoherence)")

    print("\n--- Physical Interpretation & Approximation ---")
    print("This correction is typically negative, signifying a reduction in conductivity.")
    print("\nIn the common physical limit where the phase coherence length is much larger than the elastic mean free path (L_p >> l_e), the term arctan(L_p / l_e) approaches π/2.")
    print("In this limit, the equation can be approximated as:")

    # --- Construct the approximate equation string ---
    # Δσ ≈ (e² / (2*π*ħ*L_p)) - (e² / (π²*ħ*l_e))
    
    approx_part1_numerator = f"{e}²"
    approx_part1_denominator = f"(2 * {pi} * {hbar} * {L_phi})"
    approx_part1 = f"{approx_part1_numerator} / {approx_part1_denominator}"

    approx_part2_numerator = f"{e}²"
    approx_part2_denominator = f"({pi}² * {hbar} * {l_e})"
    approx_part2 = f"{approx_part2_numerator} / {approx_part2_denominator}"

    approx_equation = f"Δσ ≈ {approx_part1} - {approx_part2}"

    print("\n" + approx_equation + "\n")
    print("This shows the correction consists of two competing parts:")
    print(f"1. A positive, temperature-dependent contribution: {approx_part1}")
    print(f"   (This term dominates the change in conductivity with temperature, as L_p is temperature-dependent)")
    print(f"2. A negative, constant offset (related to the cutoff): -{approx_part2}")

if __name__ == '__main__':
    evaluate_quantum_correction_conductivity()