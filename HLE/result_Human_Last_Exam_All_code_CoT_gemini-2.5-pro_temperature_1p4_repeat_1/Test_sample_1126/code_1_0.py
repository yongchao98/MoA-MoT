import numpy as np

# Data from the experiment
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

# --- Step 1: Calculate and print the average percentages ---
print("Analysis of Fish Larvae Behavior under different CO2 levels")
print("-" * 60)

for soundscape, conditions in data.items():
    print(f"Soundscape: {soundscape}")
    
    # Control condition
    control_data = conditions["control"]
    control_avg = np.mean(control_data)
    print(f"  Control (Current CO2):")
    print(f"    Data: {control_data}")
    print(f"    Average time spent near speaker: {control_avg:.2f}%")
    
    # High CO2 condition
    co2_data = conditions["co2"]
    co2_avg = np.mean(co2_data)
    print(f"  High CO2 (Year 2100 prediction):")
    print(f"    Data: {co2_data}")
    print(f"    Average time spent near speaker: {co2_avg:.2f}%")
    print("-" * 60)

# --- Step 2: Interpretation and Conclusion ---
print("\nConclusion based on the analysis:")

# Natural Habitat Identification
temperate_control_avg = np.mean(data["Temperate reef soundscape"]["control"])
tropical_control_avg = np.mean(data["Tropical estuarine soundscape"]["control"])

print(f"1. Natural Habitat: Under control conditions, larvae showed a strong attraction to the Tropical estuarine soundscape (average {tropical_control_avg:.2f}% > 50%), but not to the Temperate reef soundscape (average {temperate_control_avg:.2f}% < 50%). This indicates the tropical estuarine is a natural habitat.")

# Effect of High CO2
tropical_co2_avg = np.mean(data["Tropical estuarine soundscape"]["co2"])
print(f"2. Effect of High CO2 on Tropical Estuarine Settlement: Under high CO2, the attraction to the tropical estuarine soundscape was lost; the average time spent dropped from {tropical_control_avg:.2f}% to {tropical_co2_avg:.2f}%. This represents a major disturbance, meaning the larvae would not settle efficiently.")

# Evaluating the Answer Choices
print("3. Evaluating the Options:")
print("   - Option A is incorrect because settlement in the temperate reef is not 'efficient' in either scenario, and the premise about it being a natural habitat is weak.")
print("   - Option B is incorrect because the disturbance in the temperate reef is minor, and it is not the primary natural habitat.")
print("   - Option C correctly states that fish will not settle in the tropical estuarine efficiently under high CO2. This aligns perfectly with the data showing a drop from 56.33% attraction to 43.67% avoidance/neutrality.")
print("   - Option D is true but too general.")
print("   - Options E and F are incorrect because they misidentify the natural habitat or the effect of CO2.")
print("   - Option H has a weak premise about the temperate reef being a natural habitat and its conclusion is less specific than C's.")
print("\nTherefore, Option C is the best description of the experimental results.")

print("\nFinal Answer Calculation:")
print("The key comparison is for the Tropical Estuarine Soundscape.")
control_sum_str = ' + '.join(map(str, data["Tropical estuarine soundscape"]["control"]))
co2_sum_str = ' + '.join(map(str, data["Tropical estuarine soundscape"]["co2"]))

print(f"Control Average = ({control_sum_str}) / 6 = {tropical_control_avg:.2f}% (Attraction)")
print(f"High CO2 Average = ({co2_sum_str}) / 6 = {tropical_co2_avg:.2f}% (Loss of Attraction)")
print("This shows settlement in the tropical estuarine will not be as efficient.")