import numpy as np

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

# --- Analysis ---
print("Analyzing the fish larvae behavior data...")
print("-" * 40)

for soundscape, values in data.items():
    control_avg = np.mean(values["control"])
    co2_avg = np.mean(values["Co2"])

    print(f"\nSoundscape: {soundscape}")
    print(f"  - Control (Current CO2 Level):")
    print(f"    - Data: {values['control']}")
    print(f"    - Average time near speaker: {control_avg:.2f}%")
    
    # Interpretation of control behavior
    if control_avg > 55:
        print("    - Interpretation: Strong attraction (likely a natural habitat).")
    elif control_avg > 50:
        print("    - Interpretation: Mild attraction.")
    elif control_avg < 45:
        print("    - Interpretation: Repulsion.")
    else:
        print("    - Interpretation: Indifference/Neutral.")

    print(f"\n  - Elevated CO2 (Year 2100 Level):")
    print(f"    - Data: {values['Co2']}")
    print(f"    - Average time near speaker: {co2_avg:.2f}%")

    # Interpretation of CO2 behavior
    if co2_avg < 45:
         print("    - Interpretation: Repulsion/Aversion.")
    elif co2_avg > 50:
         print("    - Interpretation: Attraction/Indifference.")
    else:
        print("    - Interpretation: Indifference/Neutral.")

print("-" * 40)

# --- Conclusion based on Analysis ---
print("\nConclusion:")
print("1. Natural Habitat: Under control conditions, larvae show a strong attraction (avg 56.33%) to the 'Tropical estuarine soundscape', indicating it is a preferred natural habitat. The response to the 'Temperate reef soundscape' is neutral (avg 48.17%).")
print("2. Impact of CO2: For the 'Tropical estuarine soundscape', elevated CO2 levels cause a dramatic shift from strong attraction (56.33%) to aversion (43.67%). This means the larvae can no longer effectively navigate to their habitat.")
print("\nEvaluating the choices:")
print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
print("The second part of this statement is strongly supported by our analysis. The ability to settle in the tropical estuarine is clearly impaired. While the data doesn't strongly support the temperate reef as a primary habitat, choice C contains the most accurate conclusion based on the provided data.")

print("\nFinal Answer Choice is C.")
<<<C>>>