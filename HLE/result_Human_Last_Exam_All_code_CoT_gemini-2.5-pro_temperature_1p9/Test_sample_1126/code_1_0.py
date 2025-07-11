import numpy as np

# Store the data in a dictionary
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

print("Calculating the average time spent near the speaker for each condition:\n")

# Iterate through the data to calculate and print averages
for soundscape, conditions in data.items():
    print(f"--- {soundscape} ---")
    for condition, percentages in conditions.items():
        # Create a string of the numbers being added, e.g., "40 + 50 + ..."
        calculation_str = " + ".join(map(str, percentages))
        num_points = len(percentages)
        # Calculate the average using numpy for precision
        average = np.mean(percentages)
        # Print the full equation and the result
        print(f"Average for {condition.upper()}: ({calculation_str}) / {num_points} = {average:.2f}%")
    print("\n")

print("--- Analysis of Results ---")
print("1. Natural Habitat Identification (Control Conditions):")
print("   - The larvae show the strongest attraction to the 'Tropical estuarine soundscape' (56.33%), suggesting it is a key habitat.")
print("   - The response to the 'Temperate reef soundscape' (48.17%) is near neutral (50% = chance).")
print("\n2. Impact of High CO2:")
print("   - For the 'Tropical estuarine soundscape', the attraction significantly decreases from 56.33% to 43.67%. This indicates that larvae are less able to find this habitat.")
print("   - For the 'White noise', a sound they should ignore, their behavior flips from slight avoidance (46.17%) to attraction (53.83%), showing sensory confusion.")
print("\n3. Conclusion:")
print("   - The high CO2 level clearly disrupts the larvae's ability to distinguish and navigate to their natural habitat.")
print("   - Option C accurately states that under high CO2, the fish will not settle in the tropical estuarine environment as efficiently as before.")
print("   - This is strongly supported by the calculated drop in attraction from 56.33% to 43.67%.")

<<<C>>>