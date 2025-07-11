import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for the J=0->1 and J=1->2 rotational transitions.
    """

    # 1. Define constants
    # Given values
    r0_pm = 141.4      # bond length in pm
    k_Nm = 400.0       # force constant in N/m
    mH_amu = 1.008     # mass of Hydrogen in amu
    mBr_amu = 79.904   # mass of Bromine in amu

    # Physical constants
    amu_to_kg = 1.660539e-27  # kg/amu
    h_bar_Js = 1.0545718e-34   # reduced Planck constant in J*s
    e_C = 1.6021766e-19       # elementary charge in C (for J to eV conversion)

    # 2. Perform calculations for the centrifugal distortion constant, D
    # Convert units
    r0_m = r0_pm * 1e-12
    mH_kg = mH_amu * amu_to_kg
    mBr_kg = mBr_amu * amu_to_kg

    # Calculate reduced mass (mu)
    mu_kg = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)

    # Calculate moment of inertia (I)
    I_kg_m2 = mu_kg * r0_m**2

    # Calculate rotational constant (B) in Joules
    B_J = h_bar_Js**2 / (2 * I_kg_m2)

    # Calculate angular vibrational frequency squared (omega0^2)
    omega0_sq_s_neg2 = k_Nm / mu_kg

    # Calculate centrifugal distortion constant (D) in Joules
    # The formula is D = 4 * B^3 / (hbar^2 * omega0^2)
    D_J = (4 * B_J**3) / (h_bar_Js**2 * omega0_sq_s_neg2)

    print(f"Calculated centrifugal distortion constant, D = {D_J:.4e} J\n")

    # 3. Calculate and print the energy shifts for each transition
    
    # --- Transition 1: J = 0 to J = 1 ---
    J1 = 0
    # The magnitude of the energy shift is |ΔE_shift| = 4 * D * (J+1)^3
    shift_J1_joules = 4 * D_J * (J1 + 1)**3
    shift_J1_eV = shift_J1_joules / e_C
    shift_J1_qeV = shift_J1_eV * 1e30

    print(f"1. For the transition from J = {J1} to J = {J1+1}:")
    print(f"   The energy shift is calculated by |ΔE_shift| = 4 * D * (J+1)^3")
    print(f"   |ΔE_shift| = 4 * ({D_J:.4e} J) * ({J1}+1)^3")
    print(f"   Resulting shift = {shift_J1_joules:.4e} J")
    print(f"   In quecto-electronvolts, the shift is {shift_J1_qeV:.4f} qeV.")
    print("-" * 50)

    # --- Transition 2: J = 1 to J = 2 ---
    J2 = 1
    # The magnitude of the energy shift is |ΔE_shift| = 4 * D * (J+1)^3
    shift_J2_joules = 4 * D_J * (J2 + 1)**3
    shift_J2_eV = shift_J2_joules / e_C
    shift_J2_qeV = shift_J2_eV * 1e30

    print(f"2. For the transition from J = {J2} to J = {J2+1}:")
    print(f"   The energy shift is calculated by |ΔE_shift| = 4 * D * (J+1)^3")
    print(f"   |ΔE_shift| = 4 * ({D_J:.4e} J) * ({J2}+1)^3")
    print(f"   Resulting shift = {shift_J2_joules:.4e} J")
    print(f"   In quecto-electronvolts, the shift is {shift_J2_qeV:.4f} qeV.")
    
    # Returning the final answers for the final output tag.
    return shift_J1_qeV, shift_J2_qeV

# Run the calculation and store the results
answer1, answer2 = calculate_centrifugal_distortion_shifts()

# The final answer in the required format
# print(f"\n<<<{answer1:.4e}, {answer2:.4e}>>>")
# The above line is for developing the answer, but the final output should not have it printed.

if __name__ == '__main__':
    # This block allows the script to be executed.
    # The function call is already made above to get the answers.
    pass
