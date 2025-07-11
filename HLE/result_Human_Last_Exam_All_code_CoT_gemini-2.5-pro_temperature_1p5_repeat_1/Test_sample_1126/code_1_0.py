import collections

def analyze_fish_larvae_data():
    """
    Analyzes fish larvae behavior data to determine their natural habitat
    and the effect of elevated CO2 on their settlement.
    """
    # There is a typo in the provided data ("Day2"), which is assumed to be Day 21.
    # The percentage for "Day20" under "Temperate reef soundscape" is also clarified as CO2-51%.
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

    print("--- Data Analysis ---")
    print("A percentage > 50% indicates attraction.")
    print("A percentage < 50% indicates repulsion/avoidance.")
    print("A percentage ≈ 50% indicates indifference.\n")

    results = collections.defaultdict(dict)

    for soundscape, conditions in data.items():
        print(f"--- {soundscape} ---")
        for condition, percentages in conditions.items():
            total = sum(percentages)
            count = len(percentages)
            average = total / count
            results[soundscape][condition] = average
            
            # Create the equation string
            sum_str = " + ".join(map(str, percentages))
            equation = f"({sum_str}) / {count}"
            
            print(f"Analysis for '{condition}' condition:")
            print(f"  Equation: {equation} = {average:.2f}%")
        print("-" * (len(soundscape) + 6) + "\n")

    print("--- Analysis Summary ---")
    # Step 1: Determine natural habitat from 'control' data
    control_attraction = {s: r['control'] for s, r in results.items()}
    most_attractive_habitat = max(control_attraction, key=control_attraction.get)
    
    print(f"1. Natural Habitat Identification (Control Conditions):")
    print(f"   - Attraction to Tropical Estuarine: {results['Tropical estuarine soundscape']['control']:.2f}% (Strongest attraction)")
    print(f"   - Attraction to Temperate Reef: {results['Temperate reef soundscape']['control']:.2f}% (Slight avoidance)")
    print(f"   - Attraction to White Noise: {results['White noise']['control']:.2f}% (Slight avoidance)")
    print(f"   Conclusion: The data suggests the primary natural habitat is the 'Tropical estuarine soundscape' due to the strong attraction under normal conditions.\n")

    # Step 2: Analyze the effect of high CO2
    tropical_control_avg = results['Tropical estuarine soundscape']['control']
    tropical_co2_avg = results['Tropical estuarine soundscape']['Co2']
    
    print(f"2. Impact of Elevated CO₂ on Settlement:")
    print(f"   - Under control conditions, larvae were attracted to the Tropical Estuarine sound ({tropical_control_avg:.2f}%).")
    print(f"   - Under high CO₂ conditions, larvae avoided the Tropical Estuarine sound ({tropical_co2_avg:.2f}%).")
    print(f"   Conclusion: This reversal from attraction to avoidance indicates that elevated CO₂ will severely disrupt the ability of larvae to find and settle in their natural tropical estuarine habitat, making settlement inefficient.\n")
    
    print("--- Evaluating Answer Choices ---")
    print("Choice C states that one of the natural habitats is the tropical estuarine and that at high CO₂ levels, fish will not settle there efficiently.")
    print("This aligns perfectly with our analysis. The attraction to the tropical estuarine soundscape flips from 56.33% to 43.67%, which would lead to inefficient settlement.")

    print("\nFinal Answer Choice:")
    final_answer = "C"
    print(f"<<<{final_answer}>>>")

# Execute the analysis
analyze_fish_larvae_data()