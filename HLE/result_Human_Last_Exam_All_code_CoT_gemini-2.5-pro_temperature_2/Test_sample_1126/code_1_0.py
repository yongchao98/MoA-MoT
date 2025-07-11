import collections

def analyze_larvae_behavior():
    """
    Analyzes fish larvae behavior data to determine the effect of CO2 on habitat settlement.
    """
    # Data from the experiment
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

    print("--- Analysis of Fish Larvae Auditory Cues ---")

    for soundscape, conditions in data.items():
        print(f"\nSoundscape: {soundscape}")
        for condition, percentages in conditions.items():
            # Create the string for the equation
            numbers_str = " + ".join(map(str, percentages))
            total = sum(percentages)
            count = len(percentages)
            average = total / count
            
            # Print the full equation for calculating the average
            print(f"  Average for '{condition}': ({numbers_str}) / {count} = {average:.2f}%")

    print("\n--- Conclusion ---")
    print("Under control conditions, the larvae show a clear preference for the 'Tropical estuarine soundscape' (average > 50%), indicating it's a primary natural habitat.")
    print("Under high CO2 conditions, this preference is lost; the time spent near the 'Tropical estuarine soundscape' drops significantly below 50%.")
    print("This demonstrates that high CO2 levels will disturb settlement, making larvae unable to efficiently locate their preferred tropical habitat.")
    print("\nThis analysis supports option C, which states that the fish will not settle in the tropical estuarine efficiently under future CO2 levels.")

# Run the analysis
analyze_larvae_behavior()
print("\n<<<C>>>")
