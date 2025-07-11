def analyze_fish_settlement():
    """
    Analyzes fish larvae settlement data to determine the effect of elevated CO2.
    """
    # Data for the Tropical Estuarine Soundscape
    # Represents the percentage of time larvae spent near the speaker.
    tropical_control = [68, 63, 55, 51, 49, 52]
    tropical_co2 = [32, 37, 45, 49, 51, 48]

    # Calculate the average for the control condition (representing 2024 levels)
    control_avg = sum(tropical_control) / len(tropical_control)

    # Calculate the average for the elevated CO2 condition (representing 2100 levels)
    co2_avg = sum(tropical_co2) / len(tropical_co2)

    print("Analyzing settlement efficiency in the Tropical Estuarine soundscape:")
    print("A preference above 50% indicates efficient settlement, while a value below 50% indicates inefficient settlement or avoidance.")
    
    print("\n1. Settlement efficiency under current CO₂ levels (Control):")
    # Build the equation string for the control average calculation
    control_equation_str = f"({ ' + '.join(map(str, tropical_control)) }) / {len(tropical_control)} = {control_avg:.2f}%"
    print(f"   Calculation: {control_equation_str}")
    print("   This value indicates a preference for this habitat.")

    print("\n2. Settlement efficiency under predicted 2100 CO₂ levels:")
    # Build the equation string for the CO2 average calculation
    co2_equation_str = f"({ ' + '.join(map(str, tropical_co2)) }) / {len(tropical_co2)} = {co2_avg:.2f}%"
    print(f"   Calculation: {co2_equation_str}")
    print("   This value indicates avoidance of this habitat, meaning settlement is not efficient.")

    print("\nConclusion: The data shows that at the CO₂ level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently. This supports answer choice C.")

analyze_fish_settlement()
<<<C>>>