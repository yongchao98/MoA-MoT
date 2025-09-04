import math

def check_ph_calculation():
    """
    This function verifies the step-by-step calculation of the pH for the given mixture.
    It recalculates the final pH based on the problem's constraints and compares it
    to the provided answer.
    """
    try:
        # --- Step 1: Define initial conditions and calculate moles ---
        # The provided solution correctly identifies the initial volumes and concentrations.

        # Moles of H+ from HCl (strong acid):
        # Volume = 400 mL = 0.400 L, Concentration = 0.2 M
        moles_hcl = 0.400 * 0.2
        if not math.isclose(moles_hcl, 0.080):
            return f"Incorrect calculation of moles of H+ from HCl. Expected 0.080, but calculated {moles_hcl}."

        # Moles of CH3COOH (weak acid):
        # Volume = 500 mL = 0.500 L, Concentration = 0.1 M
        moles_ch3cooh = 0.500 * 0.1
        if not math.isclose(moles_ch3cooh, 0.050):
            return f"Incorrect calculation of moles of CH3COOH. Expected 0.050, but calculated {moles_ch3cooh}."

        # Moles of OH- from Ba(OH)2 (strong base):
        # Volume = 300 mL = 0.300 L, Concentration = 0.3 M
        # Note: Ba(OH)2 provides 2 moles of OH- per mole.
        moles_oh = 0.300 * 0.3 * 2
        if not math.isclose(moles_oh, 0.180):
            return f"Incorrect calculation of moles of OH- from Ba(OH)2. Expected 0.180, but calculated {moles_oh}."

        # --- Step 2: Perform neutralization reactions ---
        # The strong base (OH-) neutralizes all available acid protons (from both strong and weak acids).
        total_moles_acid = moles_hcl + moles_ch3cooh
        total_moles_base = moles_oh

        if not math.isclose(total_moles_acid, 0.130):
             return f"Incorrect calculation of total moles of acid. Expected 0.130, but calculated {total_moles_acid}."

        # The provided solution correctly identifies that the base is in excess.
        if total_moles_base <= total_moles_acid:
            return "Incorrect conclusion. The total moles of acid are greater than or equal to the total moles of base, but the solution claims the base is in excess."

        # Moles of OH- in excess = Total Base - Total Acid
        moles_oh_excess = total_moles_base - total_moles_acid
        if not math.isclose(moles_oh_excess, 0.050):
            return f"Incorrect calculation of excess OH- moles. Expected 0.050, but calculated {moles_oh_excess}."

        # --- Step 3: Calculate the final pH ---
        # Total volume of the solution
        total_volume = 0.500 + 0.400 + 0.300
        if not math.isclose(total_volume, 1.2):
            return f"Incorrect calculation of total volume. Expected 1.2 L, but calculated {total_volume} L."

        # Concentration of excess OH-
        conc_oh_excess = moles_oh_excess / total_volume
        if not math.isclose(conc_oh_excess, 0.04167, rel_tol=1e-4):
             return f"Incorrect calculation of [OH-]. Expected ~0.04167 M, but calculated {conc_oh_excess} M."

        # Calculate pOH
        poh = -math.log10(conc_oh_excess)
        if not math.isclose(poh, 1.38, abs_tol=0.01):
             return f"Incorrect calculation of pOH. Expected ~1.38, but calculated {poh:.2f}."

        # Calculate pH
        calculated_ph = 14.0 - poh
        
        # The final answer given is 12.62
        expected_ph = 12.62

        # Check if the calculated pH matches the expected pH within a reasonable tolerance for rounding.
        if math.isclose(calculated_ph, expected_ph, abs_tol=0.01):
            return "Correct"
        else:
            return f"Incorrect. The final calculated pH is {calculated_ph:.2f}, which does not match the provided answer of {expected_ph}."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_ph_calculation()
print(result)