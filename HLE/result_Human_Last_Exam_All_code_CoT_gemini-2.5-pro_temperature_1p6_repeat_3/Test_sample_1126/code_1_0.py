import numpy as np

# Data from the experiment
data = {
    'Temperate reef soundscape': {
        'control': [40, 50, 49, 53, 49, 48],
        'co2': [60, 50, 51, 47, 51, 52]
    },
    'Tropical estuarine soundscape': {
        'control': [68, 63, 55, 51, 49, 52],
        'co2': [32, 37, 45, 49, 51, 48]
    },
    'White noise': {
        'control': [30, 48, 47, 48, 52, 52],
        'co2': [70, 52, 53, 52, 48, 48]
    }
}

# --- Step 1 & 2: Calculate Averages ---
print("--- Analysis of Average Time Spent Near Speaker ---\n")

averages = {}
for soundscape, conditions in data.items():
    averages[soundscape] = {}
    print(f"Soundscape: {soundscape}")
    for condition, values in conditions.items():
        avg = np.mean(values)
        averages[soundscape][condition] = avg
        # Constructing the string for the equation
        equation_str = f"({ ' + '.join(map(str, values)) }) / {len(values)}"
        print(f"  Average for {condition.upper()}: {equation_str} = {avg:.2f}%")
    print("-" * 20)

# --- Step 3 & 4: Compare Conditions and Evaluate ---
print("\n--- Conclusions from Data ---\n")

# Determine natural habitat from control data
control_tropical_avg = averages['Tropical estuarine soundscape']['control']
control_temperate_avg = averages['Temperate reef soundscape']['control']

print(f"1. Determining Natural Habitat (Control Condition):")
print(f"   - Time at Tropical Estuarine: {control_tropical_avg:.2f}%")
print(f"   - Time at Temperate Reef: {control_temperate_avg:.2f}%")
print("   Conclusion: The larvae spend significantly more time near the Tropical Estuarine soundscape, indicating it is their primary natural habitat.\n")

# Determine disturbance in tropical estuarine
co2_tropical_avg = averages['Tropical estuarine soundscape']['co2']
print(f"2. Assessing CO2 Impact on Settlement in Tropical Estuarine:")
print(f"   - Time under control CO2: {control_tropical_avg:.2f}%")
print(f"   - Time under high CO2: {co2_tropical_avg:.2f}%")
print("   Conclusion: Under high CO2, the larvae's attraction to their natural habitat sound drops drastically. This means high CO2 levels disturb their ability to settle efficiently.\n")

print("--- Evaluating Answer Choices ---\n")
print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
print("\n- The first part about 'natural habitats' (plural) is plausible.")
print("- The second part, 'the fish will not settle in the tropical estuarine efficiently', is strongly supported by the data, as the preference drops from 56.33% to 43.67%. This is the most accurate and specific conclusion among the choices.")

<<<C>>>