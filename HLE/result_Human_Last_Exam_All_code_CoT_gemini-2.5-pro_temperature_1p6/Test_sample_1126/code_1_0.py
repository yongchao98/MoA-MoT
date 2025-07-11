import collections

def analyze_fish_larvae_data():
    """
    Analyzes experimental data on fish larvae behavior to determine the effect of CO2
    on their response to different soundscapes.
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

    print("--- Data Analysis ---\n")

    for soundscape, conditions in data.items():
        print(f"Analyzing: {soundscape}")
        for condition, percentages in conditions.items():
            # Create the equation string
            sum_str = " + ".join(map(str, percentages))
            # Calculate the sum and average
            total_sum = sum(percentages)
            average = total_sum / len(percentages)
            # Print the full calculation
            print(f"  - Average for '{condition}': ({sum_str}) / {len(percentages)} = {average:.2f}%")
        print("-" * 20)

    print("\n--- Interpretation of Results ---\n")
    print("1. Natural Habitat Identification (Control Conditions):")
    print("   - Tropical Estuarine: Average is 56.33%, showing a clear preference. This is very likely a natural habitat.")
    print("   - Temperate Reef: Average is 48.17%, close to 50%, indicating neutrality. It may be a less preferred or secondary habitat.")
    print("   - White Noise: Average is 46.17%, showing slight avoidance/neutrality as expected for an irrelevant sound.\n")

    print("2. Effect of Elevated CO2:")
    print("   - Tropical Estuarine: The preference drops from 56.33% (Control) to 43.67% (CO2). This is a significant shift from attraction to avoidance, showing that CO2 disrupts the ability to seek this habitat.")
    print("   - Temperate Reef: The behavior remains neutral (48.17% vs 51.83%). CO2 does not seem to significantly alter the preference for this sound.")
    print("   - White Noise: The behavior shifts from avoidance (46.17%) to preference (53.83%), particularly driven by a strong attraction (70%) on Day 16. This indicates sensory confusion.\n")
    
    print("--- Conclusion and Answer Selection ---\n")
    print("The data strongly indicates that the tropical estuarine is a natural habitat and that elevated CO2 levels will negatively impact the larvae's ability to settle there efficiently.")
    print("Option C states that both are natural habitats (plausible as they were chosen for the experiment) and that settlement in the tropical estuarine will be less efficient under high CO2. This aligns perfectly with our analysis.")
    print("\nFinal Answer Choice: C")


analyze_fish_larvae_data()
<<<C>>>