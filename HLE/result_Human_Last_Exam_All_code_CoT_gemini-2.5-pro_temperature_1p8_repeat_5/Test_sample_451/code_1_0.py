import scipy.constants as const

# Step 1: Define constants and input parameters
# Given values
r0 = 141.4e-12  # bond length in meters
k = 400.0  # force constant in N/m
mH_amu = 1.008  # mass of Hydrogen in amu
mBr_amu = 79.904  # mass of Bromine in amu

# Physical constants from scipy
amu_to_kg = const.atomic_mass
J_to_eV = const.electron_volt

# Step 2: Calculate the reduced mass (mu) in kg
mH_kg = mH_amu * amu_to_kg
mBr_kg = mBr_amu * amu_to_kg
mu = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)

# Step 3: Calculate the moment of inertia (I)
I = mu * r0**2

# Step 4: Calculate Rotational Constant (B) in Joules and angular vibrational frequency (omega)
B = const.hbar**2 / (2 * I)
omega = (k / mu)**0.5

# Step 5: Calculate the centrifugal distortion constant (D) in Joules
D = (4 * B**3) / (const.hbar * omega)**2

# Step 6 & 7: Calculate energy shifts and convert to qeV for each transition
def calculate_and_print_shift(J_initial, J_final):
    """Calculates the energy shift for a given transition and prints the details."""
    j_final_term = J_final**2 * (J_final + 1)**2
    j_initial_term = J_initial**2 * (J_initial + 1)**2
    
    prefactor = j_final_term - j_initial_term
    
    delta_E_J = -prefactor * D
    delta_E_eV = delta_E_J / J_to_eV
    delta_E_qeV = delta_E_eV * 1e30 # 1 qeV = 1e-30 eV, so 1 eV = 1e30 qeV
    
    print(f"Calculation for J={J_initial} to J={J_final} transition:")
    print("The formula for the energy shift is: ΔE = -D * [J_final² * (J_final + 1)² - J_initial² * (J_initial + 1)²]")
    print(f"For this transition, J_initial = {J_initial} and J_final = {J_final}.")
    print(f"ΔE = -D * [{J_final}²*({J_final}+1)² - {J_initial}²*({J_initial}+1)²]")
    print(f"ΔE = -D * [{j_final_term} - {j_initial_term}]")
    print(f"ΔE = -{prefactor} * D")
    print(f"The calculated centrifugal distortion constant D is {D:.4e} J.")
    print(f"So, the energy shift is ΔE = -{prefactor} * {D:.4e} J = {delta_E_J:.4e} J.")
    print(f"In quecto-electronvolts (qeV), this is: {delta_E_qeV:.4e} qeV.\n")
    return delta_E_qeV

# --- Main execution ---
# Transition 1: J = 0 to J = 1
delta_E1_qeV = calculate_and_print_shift(0, 1)

# Transition 2: J = 1 to J = 2
delta_E2_qeV = calculate_and_print_shift(1, 2)

# Final answer format
# We present the answers as a list [answer1, answer2]
final_answers = [delta_E1_qeV, delta_E2_qeV]
# <<<[-1.7659e+23, -1.4127e+24]>>>