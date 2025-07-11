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

# --- Analysis ---
print("--- Analysis of Fish Larvae Behavior ---")

# Step 1: Identify Natural Habitat by analyzing control conditions
print("\nStep 1: Analyzing behavior under control (current) conditions to identify natural habitat.")
for soundscape, values in data.items():
    control_avg = np.mean(values["control"])
    print(f"'{soundscape}': Average time spent in control conditions = {control_avg:.2f}%")

print("\nConclusion for Step 1:")
print("Under control conditions, larvae show the strongest attraction to the 'Tropical estuarine soundscape' (avg > 56%).")
print("This indicates that the tropical estuarine environment is a natural habitat.")
print("The response to the 'Temperate reef soundscape' is near 50%, indicating little to no preference.")

# Step 2: Assess the impact of elevated CO2
print("\nStep 2: Analyzing the impact of elevated CO2 on behavior.")
for soundscape, values in data.items():
    control_avg = np.mean(values["control"])
    co2_avg = np.mean(values["co2"])
    change = co2_avg - control_avg
    print(f"\nFor '{soundscape}':")
    print(f"  - Control average: {control_avg:.2f}%")
    print(f"  - High CO2 average: {co2_avg:.2f}%")
    print(f"  - Change: {change:+.2f}%")
    if soundscape == "Tropical estuarine soundscape":
        print("  - Observation: There is a significant decrease in attraction. This means the larvae's ability to find this habitat is impaired.")
    if soundscape == "White noise":
        print("  - Observation: Larvae switch from aversion/neutrality to attraction, indicating sensory confusion.")


# Step 3: Evaluate Answer Choices based on the analysis
print("\n--- Final Conclusion ---")
print("The data strongly supports that:")
print("1. The tropical estuarine environment is a natural habitat for the larvae.")
print("2. Elevated CO2 levels disrupt the larvae's sensory abilities, causing them to lose attraction to their natural habitat sound.")
print("3. This disruption means they will not be able to settle in the tropical estuarine habitat as efficiently as they do under current conditions.")
print("\nBased on this analysis, we evaluate the options. Option C states that one of the natural habitats is the tropical estuarine, and that at the higher CO2 level, fish will not settle there efficiently. This aligns perfectly with our findings.")

print("\nFinal Answer Choice: C")
<<<C>>>