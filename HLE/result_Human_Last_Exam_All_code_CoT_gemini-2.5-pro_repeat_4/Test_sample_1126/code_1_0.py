def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the effect of CO2 on habitat selection.
    """
    data = {
        "Temperate reef soundscape": {
            "control": [40, 50, 49, 53, 49, 48],
            "Co2": [60, 50, 51, 47, 51, 52]
        },
        "Tropical estuarine soundscape": {
            "control": [68, 63, 55, 51, 49, 52],
            "Co2": [32, 37, 45, 49, 51, 48]
        },
        "White noise": {
            "control": [30, 48, 47, 48, 52, 52],
            "Co2": [70, 52, 53, 52, 48, 48]
        }
    }

    print("--- Data Analysis ---")
    
    # --- Tropical Estuarine Soundscape Analysis ---
    tropical_control_data = data["Tropical estuarine soundscape"]["control"]
    tropical_co2_data = data["Tropical estuarine soundscape"]["Co2"]

    tropical_control_avg = sum(tropical_control_data) / len(tropical_control_data)
    tropical_co2_avg = sum(tropical_co2_data) / len(tropical_co2_data)

    # Building the equation strings
    tropical_control_sum_str = " + ".join(map(str, tropical_control_data))
    tropical_co2_sum_str = " + ".join(map(str, tropical_co2_data))

    print("\nAnalysis for Tropical Estuarine Soundscape:")
    print(f"Under control conditions (current CO2), the larvae were attracted to the sound.")
    print("Calculation for average time spent:")
    print(f"({tropical_control_sum_str}) / {len(tropical_control_data)} = {tropical_control_avg:.2f}%")
    
    print("\nUnder high CO2 conditions, the attraction was lost.")
    print("Calculation for average time spent:")
    print(f"({tropical_co2_sum_str}) / {len(tropical_co2_data)} = {tropical_co2_avg:.2f}%")
    
    print("\n--- Conclusion ---")
    print("The data shows that under normal conditions, the larvae are attracted to the Tropical Estuarine soundscape, suggesting it's a natural habitat.")
    print("Under high CO2 levels, this attraction is significantly reduced (from 56.33% to 43.67%).")
    print("This means the larvae will not settle in the tropical estuarine as efficiently under future CO2 conditions.")
    print("This supports answer choice C.")
    
analyze_fish_behavior()
<<<C>>>