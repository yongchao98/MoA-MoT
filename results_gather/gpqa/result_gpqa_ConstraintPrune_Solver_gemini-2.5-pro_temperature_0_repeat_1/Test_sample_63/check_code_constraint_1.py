import math

def check_chemistry_problem():
    """
    This function verifies the solution to the given chemistry problem by simulating the reactions
    and comparing the calculated results with the observed data.
    """
    # --- Given data from the problem ---
    total_mass_salts = 7.20  # g
    mass_increase_tube1_obs = 3.60  # g (H2O)
    mass_increase_tube2_obs = 0.0   # g (acidic gas)
    mass_increase_tube3_obs = 0.80  # g (O from oxidant)
    volume_gas_C_obs = 2.24  # L at STP

    # --- Constants ---
    MOLAR_MASS = {'H': 1.008, 'N': 14.007, 'O': 15.999}
    STP_MOLAR_VOLUME = 22.414  # L/mol

    # --- Proposed Answer's Salts ---
    # Salt A: Ammonium Nitrite (NH4NO2)
    # Salt B: Ammonium Nitrate (NH4NO3)
    salt_A_formula = {'N': 2, 'H': 4, 'O': 2}
    salt_B_formula = {'N': 2, 'H': 4, 'O': 3}
    proposed_total_atoms = 17

    # --- Verification Steps ---

    # 1. Calculate molar masses of the proposed salts
    molar_mass_A = sum(MOLAR_MASS[atom] * count for atom, count in salt_A_formula.items())
    molar_mass_B = sum(MOLAR_MASS[atom] * count for atom, count in salt_B_formula.items())

    # 2. Calculate moles (n) based on total mass and equimolar condition
    # n * M_A + n * M_B = total_mass_salts
    try:
        n = total_mass_salts / (molar_mass_A + molar_mass_B)
    except ZeroDivisionError:
        return "Error: Molar masses of salts are zero."

    # 3. Simulate reactions and calculate products
    # Reaction A: NH4NO2 -> N2 + 2H2O
    # Reaction B: NH4NO3 -> N2O + 2H2O
    
    # Moles of H2O produced
    moles_H2O_from_A = 2 * n
    moles_H2O_from_B = 2 * n
    total_moles_H2O = moles_H2O_from_A + moles_H2O_from_B
    
    # Moles of other gases
    moles_N2_from_A = n
    moles_N2O_from_B = n

    # 4. Check calculated results against observed data
    
    # Check 1: Mass of H2O absorbed in Tube #1
    molar_mass_H2O = 2 * MOLAR_MASS['H'] + MOLAR_MASS['O']
    mass_increase_tube1_calc = total_moles_H2O * molar_mass_H2O
    if not math.isclose(mass_increase_tube1_calc, mass_increase_tube1_obs, rel_tol=1e-2):
        return f"Incorrect: The calculated mass of H2O is {mass_increase_tube1_calc:.2f} g, but the observed mass was {mass_increase_tube1_obs} g. The constraint for tube #1 is not satisfied."

    # Check 2: Mass of acidic gas absorbed in Tube #2
    # The products are N2, N2O, and H2O. None are acidic. This condition is satisfied.

    # Check 3: Mass increase of Tube #3 due to reaction with N2O
    # Reaction: Cu + N2O -> CuO + N2. The mass increase is the mass of oxygen from N2O.
    mass_increase_tube3_calc = moles_N2O_from_B * MOLAR_MASS['O']
    if not math.isclose(mass_increase_tube3_calc, mass_increase_tube3_obs, rel_tol=1e-2):
        return f"Incorrect: The calculated mass increase of tube #3 is {mass_increase_tube3_calc:.2f} g, but the observed increase was {mass_increase_tube3_obs} g. The constraint for tube #3 is not satisfied."

    # Check 4: Volume of final gas C (N2)
    # Final N2 comes from NH4NO2 directly and from the decomposition of N2O.
    total_moles_gas_C = moles_N2_from_A + moles_N2O_from_B # 1 mole N2O produces 1 mole N2
    volume_gas_C_calc = total_moles_gas_C * STP_MOLAR_VOLUME
    if not math.isclose(volume_gas_C_calc, volume_gas_C_obs, rel_tol=1e-2):
        return f"Incorrect: The calculated volume of the final gas C is {volume_gas_C_calc:.2f} L, but the observed volume was {volume_gas_C_obs} L. The final volume constraint is not satisfied."

    # Check 5: Verify the total atom count from the proposed answer
    atoms_A = sum(salt_A_formula.values())
    atoms_B = sum(salt_B_formula.values())
    total_atoms_calc = atoms_A + atoms_B
    if total_atoms_calc != proposed_total_atoms:
        return f"Incorrect: The total number of atoms for NH4NO2 and NH4NO3 is {total_atoms_calc}, which does not match the proposed answer's value of {proposed_total_atoms}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_chemistry_problem()
print(result)