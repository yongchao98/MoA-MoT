import math
# Using scipy for precise physical constants is good practice.
# If you don't have scipy, you can define the constants manually:
# e = 1.602176634e-19  # Elementary charge in Coulombs
# hbar = 1.054571817e-34 # Reduced Planck constant in J*s
try:
    from scipy.constants import e, hbar
except ImportError:
    print("SciPy not found. Using manually defined constants.")
    e = 1.602176634e-19
    hbar = 1.054571817e-34

def calculate_quantum_correction():
    """
    Calculates the weak localization correction to conductivity in 3D.
    """
    # --- Step 1: Define material parameters ---
    # These are typical values for a doped semiconductor at low temperatures.
    # l_e: elastic mean free path (in meters)
    l_e = 10e-9  # 10 nm
    # L_phi: phase coherence length (in meters)
    # This must be greater than l_e for the effect to be observed.
    L_phi = 100e-9 # 100 nm

    print("Evaluating the quantum correction to conductivity (Δσ) in a 3D system.")
    print("Formula: Δσ = (e² / (2 * π² * ħ)) * (1/L_φ - 1/l_e)\n")

    # --- Step 2: Calculate each part of the equation ---
    prefactor = (e**2) / (2 * math.pi**2 * hbar)
    inv_L_phi = 1 / L_phi
    inv_l_e = 1 / l_e
    length_term = inv_L_phi - inv_l_e
    delta_sigma = prefactor * length_term

    # --- Step 3: Display the parameters and constants ---
    print("Constants and Parameters:")
    print(f"  e (Elementary Charge) = {e:.4e} C")
    print(f"  ħ (Reduced Planck Constant) = {hbar:.4e} J·s")
    print(f"  l_e (Elastic Mean Free Path) = {l_e:.1e} m")
    print(f"  L_φ (Phase Coherence Length) = {L_phi:.1e} m\n")

    # --- Step 4: Display the equation with numbers plugged in ---
    print("Calculation Breakdown:")
    print(f"  Prefactor (e² / (2 * π² * ħ)) = ({e:.4e}² / (2 * {math.pi:.3f}² * {hbar:.4e})) = {prefactor:.4e} S")
    print(f"  Length Term (1/L_φ - 1/l_e) = (1/{L_phi:.1e} - 1/{l_e:.1e}) = ({inv_L_phi:.2e} - {inv_l_e:.2e}) = {length_term:.2e} m⁻¹\n")

    print("Final Equation with Values:")
    print(f"Δσ = {prefactor:.4e} S * ({inv_L_phi:.2e} m⁻¹ - {inv_l_e:.2e} m⁻¹) \n")


    # --- Step 5: Print the final result ---
    print("Final Result:")
    print(f"Quantum Correction to Conductivity (Δσ) = {delta_sigma:.2f} S/m")
    print("(The negative sign indicates a reduction in conductivity, as expected for weak localization)")


if __name__ == '__main__':
    calculate_quantum_correction()
    # The final numerical answer is part of the printed output,
    # but to conform to the format request, we extract the final value.
    # Recalculating here just to capture the value for the final marker.
    l_e = 10e-9
    L_phi = 100e-9
    prefactor = (e**2) / (2 * math.pi**2 * hbar)
    length_term = (1 / L_phi) - (1 / l_e)
    final_answer = prefactor * length_term
    print(f"\n<<<{final_answer:.2f}>>>")
