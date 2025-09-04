import math

def check_chemistry_problem():
    """
    This function checks the solution to the chemistry problem by verifying all
    given constraints.
    """
    # --- 1. Define Constants and Problem Data ---
    # Molar masses (g/mol)
    MOLAR_MASS = {'H': 1.008, 'N': 14.007, 'O': 16.000}
    MOLAR_MASS_H2O = 2 * MOLAR_MASS['H'] + MOLAR_MASS['O']
    MOLAR_MASS_O2 = 2 * MOLAR_MASS['O']
    MOLAR_MASS_N2 = 2 * MOLAR_MASS['N']
    
    # Molar volume of a gas at STP (L/mol)
    MOLAR_VOLUME_STP = 22.4
    
    # Given experimental data
    initial_mass_mixture = 7.20  # g
    mass_h2o_produced = 3.60   # g
    mass_o2_produced = 0.80    # g
    volume_gas_c_produced = 2.24 # L
    
    # Tolerance for floating point comparisons
    TOLERANCE = 0.01

    # --- 2. Calculate Moles of Products from Experimental Data ---
    exp_moles_h2o = mass_h2o_produced / MOLAR_MASS_H2O
    exp_moles_o2 = mass_o2_produced / MOLAR_MASS_O2
    exp_moles_gas_c = volume_gas_c_produced / MOLAR_VOLUME_STP

    # --- 3. Verify Mass Conservation and Identify Gas C ---
    # The mass of Gas C can be found by subtracting the masses of the other products
    # from the initial mass.
    mass_gas_c = initial_mass_mixture - mass_h2o_produced - mass_o2_produced
    
    # Calculate the molar mass of Gas C to confirm its identity as N2
    molar_mass_gas_c = mass_gas_c / exp_moles_gas_c
    if not math.isclose(molar_mass_gas_c, MOLAR_MASS_N2, rel_tol=TOLERANCE):
        return (f"Incorrect: Mass conservation check failed or Gas C is not N2. "
                f"Calculated molar mass of Gas C is {molar_mass_gas_c:.2f} g/mol, "
                f"but expected N2 ({MOLAR_MASS_N2:.2f} g/mol).")
    
    exp_moles_n2 = exp_moles_gas_c

    # --- 4. Test the Hypothesis: Salts are NH4NO2 and NH4NO3 ---
    # Decomposition reactions:
    # Salt A (NH4NO2) -> N2 + 2H2O
    # Salt B (2NH4NO3) -> 2N2 + O2 + 4H2O
    
    # Since the mixture is equimolar, let 'n' be the moles of each salt.
    # Only NH4NO3 produces O2. The reaction 2NH4NO3 -> O2 means n moles of NH4NO3 produce n/2 moles of O2.
    
    # Calculate 'n' from the moles of O2
    n = exp_moles_o2 / 0.5
    
    # --- 5. Verify Stoichiometry for Other Products ---
    # Theoretical moles of N2 = n (from NH4NO2) + n (from NH4NO3) = 2n
    # Theoretical moles of H2O = 2n (from NH4NO2) + 2n (from NH4NO3) = 4n
    
    theoretical_moles_n2 = 2 * n
    theoretical_moles_h2o = 4 * n
    
    if not math.isclose(theoretical_moles_n2, exp_moles_n2, rel_tol=TOLERANCE):
        return (f"Incorrect: N2 amount is inconsistent. "
                f"Expected {theoretical_moles_n2:.3f} mol, but experimental value is {exp_moles_n2:.3f} mol.")
    
    if not math.isclose(theoretical_moles_h2o, exp_moles_h2o, rel_tol=TOLERANCE):
        return (f"Incorrect: H2O amount is inconsistent. "
                f"Expected {theoretical_moles_h2o:.3f} mol, but experimental value is {exp_moles_h2o:.3f} mol.")

    # --- 6. Verify Total Mass Constraint ---
    molar_mass_nh4no2 = 2 * MOLAR_MASS['N'] + 4 * MOLAR_MASS['H'] + 2 * MOLAR_MASS['O']
    molar_mass_nh4no3 = 2 * MOLAR_MASS['N'] + 4 * MOLAR_MASS['H'] + 3 * MOLAR_MASS['O']
    
    theoretical_total_mass = n * (molar_mass_nh4no2 + molar_mass_nh4no3)
    
    if not math.isclose(theoretical_total_mass, initial_mass_mixture, rel_tol=TOLERANCE):
        return (f"Incorrect: Total mass constraint not satisfied. "
                f"Calculated mass is {theoretical_total_mass:.2f} g, but initial mass was {initial_mass_mixture:.2f} g.")

    # --- 7. Calculate the Final Answer (Total Atoms) ---
    # Atoms in NH4NO2: 2*N + 4*H + 2*O = 8 atoms
    # Atoms in NH4NO3: 2*N + 4*H + 3*O = 9 atoms
    total_atoms = 8 + 9
    
    # The correct answer is D, which is 17.
    if total_atoms != 17:
        return f"Incorrect: The final calculation of atoms is wrong. Expected 17, but got {total_atoms}."

    return "Correct"

# Run the check
result = check_chemistry_problem()
print(result)