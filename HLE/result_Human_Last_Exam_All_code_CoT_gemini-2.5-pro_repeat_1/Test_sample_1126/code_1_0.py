import collections

def analyze_larvae_behavior():
    """
    Analyzes fish larvae behavior data to determine their natural habitat
    and the effect of elevated CO2 levels.
    """
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

    print("--- Data Analysis ---\n")

    results = collections.defaultdict(dict)

    for soundscape, conditions in data.items():
        print(f"Analyzing: {soundscape}")
        for condition, percentages in conditions.items():
            # Calculate the sum for the equation
            total_sum = sum(percentages)
            # Calculate the count for the equation
            count = len(percentages)
            # Calculate the average
            average = total_sum / count
            results[soundscape][condition] = average
            
            # Create the equation string with each number
            equation_str = " + ".join(map(str, percentages))
            
            print(f"  - {condition.upper()} Condition Average:")
            print(f"    Equation: ({equation_str}) / {count}")
            print(f"    Result: {total_sum} / {count} = {average:.2f}%\n")

    print("\n--- Conclusion from Analysis ---\n")
    
    # 1. Determine Natural Habitat
    control_averages = {s: r['control'] for s, r in results.items()}
    natural_habitat = max(control_averages, key=control_averages.get)
    print(f"1. Natural Habitat Identification:")
    print(f"   Under control conditions, the larvae showed the strongest attraction to the '{natural_habitat}' ({control_averages[natural_habitat]:.2f}%).")
    print(f"   This indicates it is their primary natural habitat.\n")
    
    # 2. Analyze CO2 effect on natural habitat
    control_attraction = results[natural_habitat]['control']
    co2_attraction = results[natural_habitat]['co2']
    print(f"2. Effect of Elevated CO2 on Natural Habitat ({natural_habitat}):")
    print(f"   The average time spent near the habitat sound dropped from {control_attraction:.2f}% (attraction) in the control group to {co2_attraction:.2f}% (avoidance) in the high CO2 group.")
    print(f"   This shows that elevated CO2 disrupts their ability to navigate to their habitat, meaning they will not settle there as efficiently.\n")

    # 3. Evaluate choices
    print("3. Evaluating Answer Choices:")
    print("   - Choice C states that one of the habitats is the tropical estuarine and that at high CO2 levels, the fish will not settle there as efficiently.")
    print("   - Our analysis confirms that the tropical estuarine is the natural habitat and that settlement efficiency is severely impacted by high CO2 (attraction drops from 56.33% to 43.67%).")
    print("   - While the claim that the temperate reef is also a natural habitat is not strongly supported by the data (control average is < 50%), the conclusion about the tropical estuarine is the most accurate and significant finding described in the options.")
    print("\nTherefore, Choice C is the best answer.")


analyze_larvae_behavior()
<<<C>>>