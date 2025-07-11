def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data under control and high CO2 conditions.
    """
    data = {
        "Temperate reef soundscape": {
            "control": [40, 50, 49, 53, 49, 48],
            "co2": [60, 50, 51, 47, 51, 52]
        },
        "Tropical estuarine soundscape": {
            "control": [68, 63, 55, 51, 49, 52],
            "co2": [32, 37, 45, 49, 51, 48]
        },
        "White noise": {
            "control": [30, 48, 47, 48, 52, 52],
            "co2": [70, 52, 53, 52, 48, 48]
        }
    }

    print("--- Analysis of Fish Larvae Auditory Response ---")

    for soundscape, conditions in data.items():
        print(f"\nSoundscape: {soundscape}")
        
        # Control condition analysis
        control_data = conditions["control"]
        control_avg = sum(control_data) / len(control_data)
        control_eq_str = " + ".join(map(str, control_data))
        print(f"Control (Current CO2):")
        print(f"  - Data: {control_data}")
        print(f"  - Average Calculation: ({control_eq_str}) / {len(control_data)} = {control_avg:.2f}%")
        if control_avg > 52:
            print("  - Interpretation: Strong attraction (likely habitat).")
        elif control_avg < 48:
            print("  - Interpretation: Avoidance (irrelevant/harmful sound).")
        else:
            print("  - Interpretation: Neutral response.")

        # High CO2 condition analysis
        co2_data = conditions["co2"]
        co2_avg = sum(co2_data) / len(co2_data)
        co2_eq_str = " + ".join(map(str, co2_data))
        print(f"High CO2 (Year 2100):")
        print(f"  - Data: {co2_data}")
        print(f"  - Average Calculation: ({co2_eq_str}) / {len(co2_data)} = {co2_avg:.2f}%")
        if co2_avg > 52:
            print("  - Interpretation: Attraction.")
        elif co2_avg < 48:
            print("  - Interpretation: Avoidance.")
        else:
            print("  - Interpretation: Neutral response.")

    print("\n--- Conclusion ---")
    print("Under normal conditions, larvae are attracted to the 'Tropical estuarine soundscape' (avg > 50%), suggesting it is their natural habitat.")
    print("Under high CO2, this attraction is lost (avg < 50%). This means they will not settle in their natural habitat as efficiently.")
    print("This matches the conclusion in option C.")

analyze_fish_behavior()
<<<C>>>