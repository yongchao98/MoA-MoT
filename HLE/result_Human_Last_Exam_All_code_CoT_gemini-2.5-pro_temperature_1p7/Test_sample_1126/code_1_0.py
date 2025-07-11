def calculate_and_print_average(name, condition, data):
    """Calculates and prints the average of a list of data points in a specific format."""
    count = len(data)
    total = sum(data)
    average = total / count
    
    # Create the string for the sum part of the equation, e.g., "68 + 63 + ..."
    sum_str = ' + '.join(map(str, data))
    
    # Print the full equation and result
    print(f"{name} ({condition}):")
    print(f"Equation: ( {sum_str} ) / {count} = {average:.2f}%")
    print("-" * 30)

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
# Calculate and print averages for all conditions
print("--- Data Analysis: Average Time Spent Near Speaker ---\n")
for soundscape, conditions in data.items():
    for condition, percentages in conditions.items():
        calculate_and_print_average(soundscape, condition, percentages)

# --- Conclusion based on Analysis ---
# 1. Natural Habitat: Under control conditions, the larvae show a clear attraction to the
#    Tropical Estuarine soundscape (average > 56%), but a neutral response to the
#    Temperate Reef (average ~48%). This indicates the Tropical Estuarine is a natural habitat.
# 2. Effect of CO2: Under elevated CO2, the attraction to the Tropical Estuarine soundscape
#    is lost and becomes avoidance (average drops from 56% to ~44%). This means the larvae
#    can no longer efficiently find their habitat.
# 3. Evaluation of Option C:
#    - "The natural habitats of the examined fish are the tropical estuarine and the temperate reef." -> This part is questionable, as only the tropical estuarine shows strong attraction.
#    - "At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024." -> This part is strongly supported by the data, as the attraction is eliminated.
#
# Comparing all options, Option C contains the most important and accurate conclusion about the
# negative impact of CO2 on the settlement behavior in the identified natural habitat, despite the
# potential inaccuracy in its first clause. It provides the best description of the key findings.
