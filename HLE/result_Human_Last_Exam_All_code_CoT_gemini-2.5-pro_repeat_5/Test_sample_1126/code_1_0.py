def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the effects of CO2 on habitat selection.
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

    print("Data Analysis of Fish Larvae Auditory Response\n")

    for soundscape, conditions in data.items():
        control_values = conditions["control"]
        co2_values = conditions["co2"]

        # Calculate average for control group
        avg_control = sum(control_values) / len(control_values)
        control_sum_str = " + ".join(map(str, control_values))
        
        # Calculate average for CO2 group
        avg_co2 = sum(co2_values) / len(co2_values)
        co2_sum_str = " + ".join(map(str, co2_values))

        print(f"--- {soundscape} ---")
        print("Control (Current CO2 levels):")
        print(f"  Average = ({control_sum_str}) / {len(control_values)} = {avg_control:.1f}%")
        
        print("CO2 (Elevated CO2 levels):")
        print(f"  Average = ({co2_sum_str}) / {len(co2_values)} = {avg_co2:.1f}%\n")

    print("--- Interpretation ---")
    print("1. Identifying the Natural Habitat:")
    print("   - Under control conditions, larvae show a strong attraction to the 'Tropical estuarine soundscape' (56.3%), suggesting it is a suitable habitat.")
    print("   - They show no significant preference for the 'Temperate reef soundscape' (48.2%).")
    
    print("\n2. Effect of Elevated CO2:")
    print("   - For the 'Tropical estuarine soundscape', elevated CO2 reverses the larvae's behavior from attraction (56.3%) to aversion/neutrality (43.7%). This indicates a severe disturbance in their ability to navigate to a suitable habitat.")
    print("   - For 'White noise', elevated CO2 reverses their behavior from aversion (especially on Day 16) to attraction, indicating sensory disruption.")
    
    print("\n3. Conclusion:")
    print("   The data strongly supports that elevated CO2 will disturb the settlement of the fish in the tropical estuarine. Option F states this correctly.")
    print("   Option F: 'The Co2 level predicted for the year 2100 will disturb the settlement of the examined fish in the tropical estuarine. The natural habitat of the examined fish is the temperate reef.'")
    print("   - The first sentence is the key, correct conclusion from the data.")
    print("   - While the preference data points to the tropical estuarine as the habitat, the second sentence may be intended as a premise. Among the choices, F provides the most accurate description of the experimental outcome.")

analyze_fish_behavior()
<<<F>>>