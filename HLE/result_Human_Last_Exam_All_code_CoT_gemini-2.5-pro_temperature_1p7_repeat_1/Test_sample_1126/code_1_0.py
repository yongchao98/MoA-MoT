def solve_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the effect of CO2 and identify the correct conclusion.
    """

    # Data from the experiment description
    # I have corrected a typo 'Day2' to 'Day21' for consistency.
    data = {
        "Temperate reef": {
            "control": [40, 50, 49, 53, 49, 48],
            "co2": [60, 50, 51, 47, 51, 52]
        },
        "Tropical estuarine": {
            "control": [68, 63, 55, 51, 49, 52],
            "co2": [32, 37, 45, 49, 51, 48]
        },
        "White noise": {
            "control": [30, 48, 47, 48, 52, 52],
            "co2": [70, 52, 53, 52, 48, 48]
        }
    }

    print("Step 1: Analyzing behavior under Control (current CO2) Conditions to determine natural habitat.")
    print("-" * 80)
    
    # Calculate average for control conditions
    control_averages = {}
    for soundscape, values in data.items():
        control_data = values["control"]
        avg = sum(control_data) / len(control_data)
        control_averages[soundscape] = avg
        # Constructing the equation string as requested
        equation_str = " + ".join(map(str, control_data))
        print(f"Average time at '{soundscape}' (Control):")
        print(f"({equation_str}) / {len(control_data)} = {avg:.2f}%")
        
    print("\nObservation from Control Data:")
    print(" - The larvae spent significantly more than 50% of their time near the 'Tropical estuarine' soundscape (56.33%). This indicates attraction and suggests it is a natural habitat.")
    print(" - The time spent near the 'Temperate reef' soundscape is close to 50% (48.17%), indicating a neutral or weakly avoided preference.")
    print(" - The time spent near 'White noise' is below 50% (46.17%), indicating avoidance, as expected for an irrelevant cue.")
    print("\nConclusion for Step 1: The primary natural habitat for these larvae is the Tropical estuarine environment.")

    print("\n\nStep 2: Analyzing the effect of high CO2 levels on behavior.")
    print("-" * 80)

    # Calculate and print CO2 averages and compare
    for soundscape, values in data.items():
        co2_data = values["co2"]
        avg_co2 = sum(co2_data) / len(co2_data)
        avg_control = control_averages[soundscape]
        
        equation_str_co2 = " + ".join(map(str, co2_data))
        print(f"Analysis for '{soundscape}':")
        print(f"  - Control Avg: {avg_control:.2f}%")
        print(f"  - High CO2 Avg: ({equation_str_co2}) / {len(co2_data)} = {avg_co2:.2f}%")
        
        if soundscape == "Tropical estuarine":
            print("  - Impact: A strong attraction (56.33%) turned into avoidance (43.67%). This is a major disturbance to their ability to find their habitat.")
        elif soundscape == "Temperate reef":
            print("  - Impact: Behavior remains neutral (around 50%).")
        else: # White noise
            print("  - Impact: A slight avoidance (46.17%) turned into an attraction (53.83%). This shows sensory confusion.")

    print("\n\nStep 3: Evaluating the answer choices based on the analysis.")
    print("-" * 80)
    print("A: Incorrect. Temperate reef is not the primary habitat, and the settlement behavior under CO2 is not the same as control.")
    print("B: Incorrect. Misidentifies the natural habitat as temperate reef.")
    print("C: Correct. Although listing the temperate reef as a natural habitat is debatable, it correctly identifies the key finding: under high CO2, the fish will not settle efficiently in the tropical estuarine (their primary habitat), as attraction drops from 56.33% to 43.67%. This is the most important conclusion from the data.")
    print("D: True, but too general. 'C' is more specific to the experimental results.")
    print("E: Incorrect. High CO2 clearly disturbs settlement in the tropical estuarine.")
    print("F: Incorrect. Misidentifies the natural habitat as temperate reef.")
    print("H: Incorrect. Similar to C, but the first clause is less specific about the disturbance.")
    
    print("\nThe most accurate and specific conclusion supported by the calculations is C.")

    print("\n<<<C>>>")

solve_fish_behavior()