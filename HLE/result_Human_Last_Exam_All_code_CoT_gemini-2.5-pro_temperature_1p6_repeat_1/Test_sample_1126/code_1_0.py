import numpy as np

def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the natural habitat
    and the effect of high CO2 levels.
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

    print("Step 1 & 2: Calculating average attraction percentages for each condition.\n")

    results = {}
    for soundscape, conditions in data.items():
        print(f"--- Analyzing {soundscape} ---")
        results[soundscape] = {}
        for condition, percentages in conditions.items():
            # I am now calculating the average for each list of percentages.
            avg = np.mean(percentages)
            results[soundscape][condition] = avg
            print(f"Condition: {condition.upper()}")
            # I am outputting each number in the equation.
            equation_str = " + ".join(map(str, percentages))
            print(f"Calculation: ({equation_str}) / {len(percentages)} = {avg:.2f}%")
        print("-" * (len(soundscape) + 16))
        print()

    print("\nStep 3: Determining the natural habitat based on 'control' conditions.")
    control_averages = {s: r['control'] for s, r in results.items()}
    habitat = max(control_averages, key=control_averages.get)
    print(f"The average attraction under control conditions is highest for '{habitat}' ({control_averages[habitat]:.2f}%).")
    print("An average attraction greater than 50% indicates a preference for that soundscape.")
    print(f"Therefore, the 'Tropical estuarine soundscape' is the natural habitat.\n")
    
    print("Step 4: Analyzing the effect of high CO2 on settlement in the natural habitat.")
    control_habitat_avg = results[habitat]['control']
    co2_habitat_avg = results[habitat]['co2']
    print(f"In the '{habitat}', the average attraction changed from {control_habitat_avg:.2f}% (control) to {co2_habitat_avg:.2f}% (CO2).")
    print("This shows a significant decrease in attraction, turning preference into avoidance. This means the fish will not settle as efficiently.\n")

    print("Step 5: Evaluating the answer choices.")
    print("Based on the analysis, we need to find the choice that best matches these findings.")
    print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
    print("The second part of this statement, '...the fish will not settle in the tropical estuarine efficiently as in the year 2024,' is strongly supported by our analysis showing a drop from 56.33% to 43.67%.")
    print("Although the first part of the statement incorrectly includes 'temperate reef' as a habitat (our data shows it's neutral), this option captures the most critical finding of the experiment. Compared to other options which contain more significant errors, C is the best fit.")

analyze_fish_behavior()
<<<C>>>