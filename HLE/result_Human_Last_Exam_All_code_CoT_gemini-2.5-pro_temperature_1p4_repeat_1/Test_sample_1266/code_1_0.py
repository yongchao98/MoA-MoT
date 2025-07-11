def simulate_aldh_change():
    """
    This script simulates the change in ALDH levels in RAW 264.7 cells
    when treated with electrophilic compounds, based on established
    biological principles.
    """

    # --- Initial Parameters ---
    # Baseline level of ALDH protein in untreated cells (arbitrary units)
    baseline_aldh = 100.0
    # Concentration of compounds used for treatment
    concentration_uM = 50

    # --- Biological Factors (based on the rationale) ---
    # (2E)-4-Hydroxy-2-nonen-8-ynal (HNY) increases ALDH levels. We'll model this as a fold-increase.
    hny_increase_factor = 2.5
    # 4-octyl itaconate (4-OI) is a more potent activator, so its factor is higher.
    four_oi_increase_factor = 3.5

    # --- Calculations ---
    # Calculate the new ALDH amount after HNY treatment
    final_aldh_hny = baseline_aldh * hny_increase_factor
    # Calculate the new ALDH amount after 4-OI treatment
    final_aldh_4oi = baseline_aldh * four_oi_increase_factor

    # --- Output Results ---
    print("Simulation of ALDH change in RAW 264.7 cells.")
    print("-" * 50)
    print(f"Initial ALDH level: {baseline_aldh} units.")
    print("The protein responsible for sensing these compounds is Keap1.")
    print("-" * 50)

    print(f"Treatment with {concentration_uM} uM (2E)-4-Hydroxy-2-nonen-8-ynal:")
    # Outputting the equation with numbers
    print(f"Equation: {baseline_aldh} * {hny_increase_factor} = {final_aldh_hny}")
    print(f"The amount of ALDH will INCREASE to {final_aldh_hny} units.")

    print("-" * 50)
    print(f"Treatment with {concentration_uM} uM 4-OI:")
    # Outputting the equation with numbers
    print(f"Equation: {baseline_aldh} * {four_oi_increase_factor} = {final_aldh_4oi}")
    print(f"The change is MORE, increasing ALDH to {final_aldh_4oi} units.")
    print("-" * 50)

    print("\nSummary of findings:")
    print("1. ALDH Amount: Increases")
    print("2. Change with 4-OI vs HNY: More")
    print("3. Protein Involved: Keap1")


# Execute the simulation
simulate_aldh_change()