import numpy as np

# Data from the problem description
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

print("--- Analysis of Fish Larvae Behavior ---")

# 1. Identify natural habitat from control conditions
print("\nStep 1: Identifying Natural Habitat (Control Conditions)")
for soundscape, values in data.items():
    control_avg = np.mean(values["control"])
    print(f"Average time near '{soundscape}' (Control): {control_avg:.2f}%")

print("\nConclusion for Step 1:")
print("Under control conditions, the larvae show the strongest attraction to the 'Tropical estuarine soundscape' (average > 50%).")
print("They show indifference or slight avoidance to the 'Temperate reef soundscape' (average < 50%). Therefore, the tropical estuarine is the most likely primary natural habitat.")


# 2. Assess the impact of high CO2
print("\nStep 2: Assessing the Impact of High CO2")
for soundscape, values in data.items():
    control_avg = np.mean(values["control"])
    co2_avg = np.mean(values["Co2"])
    change = co2_avg - control_avg
    print(f"\nSoundscape: {soundscape}")
    print(f"  - Control Avg: {control_avg:.2f}%")
    print(f"  - High CO2 Avg: {co2_avg:.2f}%")
    print(f"  - Change in preference: {change:+.2f}%")

print("\nConclusion for Step 2:")
print("For the 'Tropical estuarine soundscape', the larvae's behavior flips from strong attraction (56.33%) to avoidance (43.67%).")
print("This indicates that at the CO2 level predicted for 2100, the larvae will not be able to effectively use this sound cue to find and settle in their habitat.")
print("Therefore, their settlement in the tropical estuarine will not be as efficient as it is today.")

# 3. Evaluate the options
print("\nStep 3: Evaluating Answer Choices")
print("Based on our analysis:")
print(" - The 'Tropical estuarine' is the clear natural habitat.")
print(" - High CO2 levels cause attraction to this habitat to turn into avoidance.")
print(" - This means settlement in the tropical estuarine will be much less efficient.")
print("\nChoice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
print("The second part of this statement is strongly supported by the data and represents the most significant finding. The first part, while debatable for the temperate reef, is plausible if we consider it a secondary or potential habitat. Among the given options, this is the most accurate and comprehensive conclusion.")

print("\n--- Final Calculation Summary ---")
# Tropical Estuarine Calculation for the final answer
tropical_control = data["Tropical estuarine soundscape"]["control"]
tropical_co2 = data["Tropical estuarine soundscape"]["Co2"]
tropical_control_avg = np.mean(tropical_control)
tropical_co2_avg = np.mean(tropical_co2)

print("For the Tropical Estuarine Soundscape:")
print(f"The average time spent near the speaker under control conditions was ({' + '.join(map(str, tropical_control))}) / {len(tropical_control)} = {tropical_control_avg:.2f}%.")
print(f"The average time spent near the speaker under high CO2 conditions was ({' + '.join(map(str, tropical_co2))}) / {len(tropical_co2)} = {tropical_co2_avg:.2f}%.")
print("This shift from attraction (>50%) to avoidance (<50%) demonstrates that settlement will not be efficient.")