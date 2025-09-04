import math

def check_titration_answer():
    """
    Checks the correctness of the given answer for the titration problem.
    """
    # --- Given constants and initial values ---
    Ka = 1.85e-5
    Kw = 1.0e-14
    initial_vol_hac = 0.020  # L (20.00 cm3)
    initial_conc_hac = 0.05   # M
    vol_water_added = 0.020  # L (20.00 cm3)
    conc_naoh = 0.1          # M

    # The answer to check (Option A)
    expected_ph_25 = 4.26
    expected_ph_equiv = 8.52

    # --- Part 1: Calculate pH at 25% titration ---
    pKa = -math.log10(Ka)
    # At 25% titration, [A-]/[HA] = 25/75 = 1/3
    ratio_A_HA = 1/3
    calculated_ph_25 = pKa + math.log10(ratio_A_HA)

    # --- Part 2: Calculate pH at the equivalence point ---
    # Initial moles of acetic acid
    moles_hac = initial_vol_hac * initial_conc_hac

    # Volume of NaOH needed to reach equivalence
    vol_naoh_equiv = moles_hac / conc_naoh

    # Total volume at equivalence point
    total_vol_equiv = initial_vol_hac + vol_water_added + vol_naoh_equiv

    # Concentration of acetate [Ac-] at equivalence
    conc_ac_equiv = moles_hac / total_vol_equiv

    # Calculate Kb for acetate
    Kb = Kw / Ka

    # Calculate [OH-] from hydrolysis: Kb = [OH-]^2 / [Ac-]
    # [OH-] = sqrt(Kb * [Ac-])
    conc_oh = math.sqrt(Kb * conc_ac_equiv)

    # Calculate pOH and then pH
    pOH = -math.log10(conc_oh)
    calculated_ph_equiv = 14.0 - pOH

    # --- Part 3: Compare calculated values with the expected answer ---
    # Use a small tolerance for floating-point comparisons
    tolerance = 0.01

    # Check pH at 25% titration
    if not math.isclose(calculated_ph_25, expected_ph_25, abs_tol=tolerance):
        return (f"Incorrect. The pH at 25% titration is incorrect. "
                f"Calculated value is {calculated_ph_25:.2f}, but the answer provides {expected_ph_25}.")

    # Check pH at equivalence point
    if not math.isclose(calculated_ph_equiv, expected_ph_equiv, abs_tol=tolerance):
        return (f"Incorrect. The pH at the equivalence point is incorrect. "
                f"Calculated value is {calculated_ph_equiv:.2f}, but the answer provides {expected_ph_equiv}.")

    # If both checks pass
    return "Correct"

# Run the check and print the result
result = check_titration_answer()
print(result)