def model_aldh_change():
    """
    This function models the relative change in ALDH levels in response to two
    different Nrf2 activators, HNY and 4-OI, based on their known potencies.
    """
    # 1. Define initial and hypothetical parameters
    baseline_aldh_units = 100.0
    # HNY is a known Nrf2 activator. Let's assign it a hypothetical potency factor.
    hny_potency_factor = 2.5
    # 4-OI is known to be a more potent Nrf2 activator than HNY.
    four_oi_potency_factor = 4.0

    print("--- Modeling ALDH Change Based on Activator Potency ---")
    print(f"Baseline ALDH level: {baseline_aldh_units} units")
    print(f"Assumed Potency Factor for HNY: {hny_potency_factor}")
    print(f"Assumed Potency Factor for 4-OI: {four_oi_potency_factor} (higher because it's more potent)")
    print("-" * 55)

    # 2. Calculate the change for each compound
    # Equation: Change = (Baseline * Factor) - Baseline
    change_with_hny = (baseline_aldh_units * hny_potency_factor) - baseline_aldh_units
    change_with_4oi = (baseline_aldh_units * four_oi_potency_factor) - baseline_aldh_units

    # 3. Print the results and the final conclusion
    print("Step 1: Both HNY and 4-OI activate the Keap1-Nrf2 pathway, causing an INCREASE in ALDH.")
    print("\nCalculating the change with HNY:")
    print(f"Final Amount = {baseline_aldh_units} * {hny_potency_factor} = {baseline_aldh_units * hny_potency_factor}")
    print(f"Equation for change: ({baseline_aldh_units} * {hny_potency_factor}) - {baseline_aldh_units} = {change_with_hny} units increase")


    print("\nCalculating the change with 4-OI:")
    print(f"Final Amount = {baseline_aldh_units} * {four_oi_potency_factor} = {baseline_aldh_units * four_oi_potency_factor}")
    print(f"Equation for change: ({baseline_aldh_units} * {four_oi_potency_factor}) - {baseline_aldh_units} = {change_with_4oi} units increase")

    print("\n--- Final Conclusion ---")
    comparison = "more" if change_with_4oi > change_with_hny else "less"
    print(f"The ALDH change is an 'increase'.")
    print(f"The change with 4-OI ({change_with_4oi}) is '{comparison}' than with HNY ({change_with_hny}).")
    print("The protein involved in sensing these compounds is 'Keap1'.")
    print("\nThis corresponds to: increase, more, Keap1.")

# Execute the modeling function
model_aldh_change()