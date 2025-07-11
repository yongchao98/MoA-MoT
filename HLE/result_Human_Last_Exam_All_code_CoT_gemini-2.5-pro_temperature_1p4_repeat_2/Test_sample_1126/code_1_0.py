import collections

# Data from the experiment
data = {
    'Temperate reef soundscape': {
        'control': [40, 50, 49, 53, 49, 48],
        'Co2': [60, 50, 51, 47, 51, 52]
    },
    'Tropical estuarine soundscape': {
        'control': [68, 63, 55, 51, 49, 52],
        'Co2': [32, 37, 45, 49, 51, 48]
    },
    'White noise': {
        'control': [30, 48, 47, 48, 52, 52],
        'Co2': [70, 52, 53, 52, 48, 48]
    }
}

print("--- Data Analysis ---\n")

# Calculate and print the average for each category
for soundscape, conditions in data.items():
    print(f"Soundscape: {soundscape}")
    for condition, percentages in conditions.items():
        total = sum(percentages)
        count = len(percentages)
        average = total / count
        
        # Format the equation string
        equation_str = " + ".join(map(str, percentages))
        
        print(f"  Condition: {condition.upper()}")
        print(f"  Calculation: ({equation_str}) / {count} = {average:.2f}%")
    print("-" * 30)

print("\n--- Conclusion ---\n")
print("1. Identifying Natural Habitat: Under 'control' conditions, the average time spent near the 'Tropical estuarine soundscape' speaker is ~56.33%, which is significantly higher than for the other soundscapes. This strong attraction indicates it is a primary natural habitat for the larvae.")
print("\n2. Impact of Elevated CO2: For the 'Tropical estuarine soundscape', the average time spent near the speaker drops dramatically from 56.33% (control) to 43.67% (CO2). This shows that the larvae's attraction to their natural habitat sound is significantly impaired, meaning they will not settle as efficiently.")
print("\n3. Evaluating Options: Option C states that at elevated CO2, the fish will not settle in the tropical estuarine efficiently. This aligns perfectly with our analysis. While the claim about the temperate reef also being a natural habitat is less strongly supported by the data (average is near 50%), the core conclusion of option C is the most accurate and significant finding from the experiment.")

print("\n<<<C>>>")