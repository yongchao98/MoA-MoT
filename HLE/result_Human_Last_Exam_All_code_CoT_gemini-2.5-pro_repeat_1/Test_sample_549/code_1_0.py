import math

def evaluate_quantum_correction():
    """
    This function prints the step-by-step evaluation of the quantum correction
    to conductivity in a 3D bulk semiconductor.
    """

    # Define the symbols used in the equations for clarity
    delta_sigma = "δσ"
    e = "e"
    h_bar = "ħ"
    pi = "π"
    D = "D"
    q = "q"
    l = "l"
    tau_phi = "τ_φ"
    L_phi = "L_φ"

    # Step 1: Explain the physical origin
    print("Step 1: Physical Origin - Weak Localization")
    print("The quantum correction to conductivity in a disordered bulk semiconductor arises from the weak localization effect.")
    print("This is a quantum interference phenomenon where an electron moving on a closed-loop path interferes")
    print("constructively with its time-reversed counterpart. This enhances the probability of the electron returning")
    print("to its starting point, effectively 'localizing' it slightly and reducing the overall conductivity.")
    print("-" * 70)

    # Step 2: Mathematical Formulation
    print("Step 2: Mathematical Formulation")
    print("The correction is calculated by integrating the 'Cooperon' propagator, C(q), which represents the sum")
    print("of probabilities for all pairs of time-reversed paths. The formula relating the conductivity")
    print("correction (δσ) to the Cooperon is:")
    print(f"  {delta_sigma} = - (2 * {e}² / {h_bar}) * {D} * Integral[ d³{q} / (2{pi})³ ] * C(q)")
    print(f"In 3D, the Cooperon is C(q) = 1 / ({D}{q}² + 1/{tau_phi}), where:")
    print(f"  {D}   is the electron diffusion constant.")
    print(f"  {tau_phi} is the phase coherence time (the time over which interference is maintained).")
    print("-" * 70)

    # Step 3: Evaluating the Integral
    print("Step 3: Evaluating the Integral")
    print("We integrate over momentum q in 3D spherical coordinates. The diffusion model is valid for momentum")
    print(f"transfers q < 1/l, where l is the elastic mean free path. So, we integrate q from 0 to 1/l.")
    print(f"  Integral = Integral from 0 to 1/{l} of [ (4{pi}{q}²) / (8{pi}³) * (1 / ({D}{q}² + 1/{tau_phi})) ] d{q}")
    print(f"This simplifies to:")
    print(f"  Integral = (1 / (2{pi}²)) * Integral from 0 to 1/{l} of [ {q}² / ({D}{q}² + 1/{tau_phi}) ] d{q}")
    print("\nSolving this definite integral gives:")
    print(f"  Integral = (1 / (2{pi}²{D})) * [ 1/{l} - (1/{L_phi}) * arctan({L_phi}/{l}) ]")
    print(f"Here, {L_phi} = sqrt({D}*{tau_phi}) is the phase coherence length.")
    print("-" * 70)

    # Step 4: The Final Equation
    print("Step 4: Final Result for Quantum Correction")
    print("Substituting the result of the integral back into the formula for δσ:")
    print(f"  {delta_sigma} = - (2 * {e}² * {D} / {h_bar}) * (1 / (2{pi}²{D})) * [ 1/{l} - (1/{L_phi}) * arctan({L_phi}/{l}) ]")
    print("\nThe diffusion constants D cancel out, giving the final equation:")
    print("\n" + "="*25 + " FINAL EQUATION " + "="*25)
    final_equation = f"  {delta_sigma} = -({e}² / ({pi}² * {h_bar})) * [1/{l} - (1/{L_phi}) * arctan({L_phi} / {l})]"
    print(final_equation)
    print("="*68 + "\n")
    print("This formula shows the reduction in conductivity (note the negative sign). It depends on fundamental")
    print(f"constants ({e}, {h_bar}, {pi}) and material-specific length scales ({l} and {L_phi}).")

if __name__ == '__main__':
    evaluate_quantum_correction()
