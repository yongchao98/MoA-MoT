import math

def check_titration_ph():
    """
    This function checks the correctness of the calculated pH values for a titration problem.
    It recalculates the pH at 25% titration and at the equivalence point based on the problem statement
    and compares them to the values given in the selected answer.
    """
    try:
        # --- Problem Parameters ---
        initial_acid_conc = 0.05  # M
        initial_acid_vol = 20.00 / 1000  # L
        water_vol = 20.00 / 1000  # L
        base_conc = 0.1  # M (NaOH)
        Ka_acid = 1.85e-5
        Kw = 1.0e-14

        # --- Answer to Check ---
        # The provided answer is 'A', which corresponds to the values (4.26, 8.52).
        expected_ph_25_percent = 4.26
        expected_ph_equiv_point = 8.52

        # --- Calculation Step 1: Initial Dilution ---
        # Calculate the concentration of acetic acid after dilution.
        diluted_acid_vol = initial_acid_vol + water_vol
        initial_moles_acid = initial_acid_conc * initial_acid_vol
        diluted_acid_conc = initial_moles_acid / diluted_acid_vol

        # Check if the calculated diluted concentration matches the one in the explanation (0.025 M)
        if not math.isclose(diluted_acid_conc, 0.025, rel_tol=1e-4):
            return f"Incorrect: The calculation of the diluted acetic acid concentration is wrong. Expected 0.025 M, but calculated {diluted_acid_conc:.4f} M."

        # --- Calculation Step 2: pH at 25% Titration ---
        # Use the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
        pKa = -math.log10(Ka_acid)
        
        # At 25% titration, the ratio of [A-]/[HA] is 25/75 or 1/3.
        ratio_base_acid = 25.0 / 75.0
        calculated_ph_25_percent = pKa + math.log10(ratio_base_acid)

        # --- Calculation Step 3: pH at the Equivalence Point ---
        # Moles of base needed = initial moles of acid
        vol_base_added = initial_moles_acid / base_conc
        
        # Total volume at the equivalence point
        total_vol_equiv = diluted_acid_vol + vol_base_added
        
        # Concentration of the conjugate base (acetate) at the equivalence point
        acetate_conc = initial_moles_acid / total_vol_equiv

        # Check if the calculated acetate concentration matches the one in the explanation (0.02 M)
        if not math.isclose(acetate_conc, 0.02, rel_tol=1e-4):
            return f"Incorrect: The calculation of the acetate concentration at the equivalence point is wrong. Expected 0.02 M, but calculated {acetate_conc:.4f} M."

        # Calculate Kb for the conjugate base
        Kb_acetate = Kw / Ka_acid
        
        # Calculate [OH-] from hydrolysis: Kb = [OH-]^2 / [A-]
        hydroxide_conc = math.sqrt(Kb_acetate * acetate_conc)
        
        # Calculate pOH and then pH
        pOH = -math.log10(hydroxide_conc)
        calculated_ph_equiv_point = 14.0 - pOH

        # --- Final Verification ---
        # Compare the calculated pH values with the expected values from the answer.
        # We round to two decimal places as is standard for this type of problem.
        
        rounded_ph_25 = round(calculated_ph_25_percent, 2)
        rounded_ph_equiv = round(calculated_ph_equiv_point, 2)

        if not math.isclose(rounded_ph_25, expected_ph_25_percent, abs_tol=0.01):
            return f"Incorrect: The calculated pH at 25% titration is {rounded_ph_25}, but the answer states it should be {expected_ph_25_percent}."

        if not math.isclose(rounded_ph_equiv, expected_ph_equiv_point, abs_tol=0.01):
            return f"Incorrect: The calculated pH at the equivalence point is {rounded_ph_equiv}, but the answer states it should be {expected_ph_equiv_point}."

        return "Correct"

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Execute the check
result = check_titration_ph()
print(result)