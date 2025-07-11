import math

def evaluate_quantum_correction():
    """
    Calculates the quantum correction to conductivity for an electron in a
    3D disordered medium (bulk semiconductor) based on the weak localization theory.
    """
    print("This script evaluates the quantum correction to conductivity (δσ) in a 3D semiconductor.")
    print("The formula from weak localization theory is:")
    print("δσ = - (e² / (2 * π² * ħ)) * (1/l - 1/l_φ)\n")

    # --- Physical Constants ---
    e = 1.60217663e-19      # Elementary charge in Coulombs
    hbar = 1.05457182e-34   # Reduced Planck constant in J·s
    pi = math.pi

    # --- Assumed Parameters for a Doped Semiconductor at Low Temperature ---
    # These values can vary significantly depending on the material, doping, and temperature.
    l = 50e-9      # Elastic mean free path in meters (50 nm)
    l_phi = 500e-9 # Phase coherence length in meters (500 nm)

    print("--- Values Used in the Calculation ---")
    print(f"e (elementary charge)      = {e:.4e} C")
    print(f"ħ (reduced Planck constant)  = {hbar:.4e} J·s")
    print(f"π (pi)                       = {pi:.5f}")
    print(f"l (elastic mean free path)   = {l:.1e} m")
    print(f"l_φ (phase coherence length) = {l_phi:.1e} m\n")

    # --- Step-by-step Calculation ---
    # 1. Calculate the universal prefactor, which has units of conductance (Siemens)
    prefactor = e**2 / (2 * pi**2 * hbar)

    # 2. Calculate the length-dependent term, which has units of inverse meters (m⁻¹)
    length_term = (1/l - 1/l_phi)

    # 3. Calculate the final correction to conductivity
    delta_sigma = -prefactor * length_term

    # --- Final Equation with Numbers ---
    print("--- Final Equation ---")
    print("Plugging the values into the formula:")
    print("δσ = - (e² / (2 * π² * ħ)) * (1/l - 1/l_φ)")
    print(f"δσ = - ({prefactor:.5e} S) * ( (1 / {l:.1e} m) - (1 / {l_phi:.1e} m) )")
    print(f"δσ = - ({prefactor:.5e} S) * ({length_term:.2e} m⁻¹)\n")


    # --- Final Result ---
    print("--- Result ---")
    print(f"The calculated quantum correction to conductivity is: δσ = {delta_sigma:.2f} S/m")
    print("(This negative value indicates a reduction in conductivity, as expected from weak localization).")

    return delta_sigma

# Run the evaluation and store the final answer
final_answer = evaluate_quantum_correction()
# The final answer is wrapped in <<<>>> as requested.
# The numeric value is returned for this purpose.
print(f"\n<<<{final_answer:.2f}>>>")
