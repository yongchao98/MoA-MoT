def calculate_aldh_change():
    """
    This function models the change in ALDH levels in response to two different
    Nrf2 activators, (2E)-4-Hydroxy-2-nonen-8-ynal (HNEY) and 4-octyl itaconate (4-OI).
    """
    # Let's assume a baseline level of ALDH in the cells.
    baseline_aldh_units = 100.0

    # Let's assign hypothetical fold-increase factors for each compound.
    # 4-OI is a more potent activator, so it gets a higher factor.
    hney_increase_factor = 2.5
    four_oi_increase_factor = 4.0

    # Calculate the new ALDH levels after treatment.
    final_aldh_with_hney = baseline_aldh_units * hney_increase_factor
    final_aldh_with_4oi = baseline_aldh_units * four_oi_increase_factor

    print("Modeling ALDH Change in RAW 264.7 Cells:")
    print("-" * 40)
    print(f"Baseline ALDH level: {baseline_aldh_units} units")
    print("\nTreatment with 50 uM (2E)-4-Hydroxy-2-nonen-8-ynal (HNEY):")
    # Outputting each number in the final equation as requested.
    print(f"Calculation: {baseline_aldh_units} * {hney_increase_factor} = {final_aldh_with_hney} units")

    print("\nTreatment with 50 uM 4-OI:")
    # Outputting each number in the final equation as requested.
    print(f"Calculation: {baseline_aldh_units} * {four_oi_increase_factor} = {final_aldh_with_4oi} units")
    print("-" * 40)

    # State the final conclusion based on the model.
    print("\nConclusion:")
    print("1. ALDH levels INCREASE in response to both electrophiles.")
    print(f"2. The increase with 4-OI ({final_aldh_with_4oi} units) is MORE than with HNEY ({final_aldh_with_hney} units).")
    print("3. This signaling process is mediated by the KEAP1-Nrf2 pathway.")

# Run the calculation and print the results.
calculate_aldh_change()