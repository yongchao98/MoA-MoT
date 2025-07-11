import collections

def analyze_fish_larvae_data():
    """
    Analyzes experimental data on fish larvae's response to soundscapes
    under different CO2 conditions.
    """
    # There is a typo in the original problem description "Day2".
    # Based on the sequence, it is corrected to "Day 21".
    data = {
        "Temperate reef soundscape": {
            "control": collections.OrderedDict([("Day 16", 40), ("Day 17", 50), ("Day 18", 49), ("Day 19", 53), ("Day 20", 49), ("Day 21", 48)]),
            "Co2": collections.OrderedDict([("Day 16", 60), ("Day 17", 50), ("Day 18", 51), ("Day 19", 47), ("Day 20", 51), ("Day 21", 52)])
        },
        "Tropical estuarine soundscape": {
            "control": collections.OrderedDict([("Day 16", 68), ("Day 17", 63), ("Day 18", 55), ("Day 19", 51), ("Day 20", 49), ("Day 21", 52)]),
            "Co2": collections.OrderedDict([("Day 16", 32), ("Day 17", 37), ("Day 18", 45), ("Day 19", 49), ("Day 20", 51), ("Day 21", 48)])
        },
        "White noise": {
            "control": collections.OrderedDict([("Day 16", 30), ("Day 17", 48), ("Day 18", 47), ("Day 19", 48), ("Day 20", 52), ("Day 21", 52)]),
            "Co2": collections.OrderedDict([("Day 16", 70), ("Day 17", 52), ("Day 18", 53), ("Day 19", 52), ("Day 20", 48), ("Day 21", 48)])
        }
    }

    print("Analyzing Experimental Data:\n")
    
    results = {}

    for soundscape, conditions in data.items():
        print(f"--- {soundscape} ---")
        results[soundscape] = {}
        for condition, daily_data in conditions.items():
            values = list(daily_data.values())
            average = sum(values) / len(values)
            results[soundscape][condition] = average
            
            # Create the string for the equation, showing each number
            values_str = ' + '.join(map(str, values))
            print(f"Average for '{condition}': ({values_str}) / {len(values)} = {average:.2f}%")
        
        print("-" * (len(soundscape) + 6) + "\n")

    print("--- Conclusions from Analysis ---\n")
    
    # 1. Determine natural habitat from control data
    control_tropical_avg = results['Tropical estuarine soundscape']['control']
    control_temperate_avg = results['Temperate reef soundscape']['control']
    
    print(f"1. Natural Habitat Inference (Control Conditions):")
    print(f"   - Attraction to Tropical Estuarine Sound: {control_tropical_avg:.2f}%")
    print(f"   - Attraction to Temperate Reef Sound: {control_temperate_avg:.2f}%")
    print("   The significantly higher attraction to the tropical estuarine soundscape suggests it is a primary natural habitat. The temperate reef may also be a suitable habitat.")
    
    # 2. Analyze the effect of elevated CO2
    co2_tropical_avg = results['Tropical estuarine soundscape']['Co2']
    change_tropical = co2_tropical_avg - control_tropical_avg

    print(f"\n2. Effect of Elevated CO2:")
    print(f"   - In the 'Tropical estuarine soundscape', the average time spent near the speaker dropped from {control_tropical_avg:.2f}% to {co2_tropical_avg:.2f}%.")
    print("     This shows that elevated CO2 levels significantly impair the larvae's ability to navigate towards their primary habitat sound.")
    
    print("\n--- Final Evaluation ---\n")
    print("Based on the analysis, the elevated CO2 disrupts the larvae's sensory ability, causing them to no longer be attracted to their primary habitat sound (tropical estuarine).")
    print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
    print("This statement accurately reflects our findings.")

analyze_fish_larvae_data()
print("<<<C>>>")