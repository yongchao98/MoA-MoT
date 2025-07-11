def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the effect of CO2.
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

    print("--- Calculating Average Time Spent Near Speaker ---\n")

    for soundscape, conditions in data.items():
        for condition, values in conditions.items():
            total = sum(values)
            count = len(values)
            average = total / count
            # Create the sum string like "40 + 50 + ..."
            sum_str = ' + '.join(map(str, values))
            print(f"Average for {soundscape} ({condition}):")
            print(f"({sum_str}) / {count} = {average:.2f}%")
            print("-" * 40)
    
    print("\n--- Analysis Summary ---")
    print("1. Under control conditions, larvae showed a clear preference for the 'Tropical estuarine soundscape' (56.33%) over the 'Temperate reef soundscape' (48.17%), suggesting the former is their primary habitat.")
    print("2. Under elevated CO2, the preference for the 'Tropical estuarine soundscape' was lost, dropping to 43.67% (an aversion).")
    print("3. This indicates that at the CO2 level predicted for 2100, the fish will not settle in the tropical estuarine efficiently.")
    print("\nBased on this analysis, the most accurate answer is C.")

analyze_fish_behavior()
<<<C>>>