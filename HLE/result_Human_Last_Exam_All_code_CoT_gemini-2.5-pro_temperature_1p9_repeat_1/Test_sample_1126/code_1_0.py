def solve_fish_behavior_problem():
    """
    Analyzes fish larvae behavior data to determine the impact of CO2
    and choose the best-fitting conclusion from a list of options.
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

    print("Step 1: Calculating average time spent near the speaker for each condition.\n")
    
    averages = {}
    for soundscape, conditions in data.items():
        print(f"--- {soundscape} ---")
        averages[soundscape] = {}
        for condition, percentages in conditions.items():
            avg = sum(percentages) / len(percentages)
            averages[soundscape][condition] = avg
            # Building the equation string with each number
            equation_str = f"({ ' + '.join(map(str, percentages)) }) / {len(percentages)}"
            print(f"'{condition.capitalize()}' Average: {equation_str} = {avg:.2f}%")
        print()

    print("Step 2: Analyzing the results.\n")

    # Analysis of natural habitat
    print("Analysis of Natural Habitat (Control Condition):")
    print(f" - Temperate reef: {averages['Temperate reef soundscape']['control']:.2f}%. This is close to 50%, suggesting indifference or slight aversion.")
    print(f" - Tropical estuarine: {averages['Tropical estuarine soundscape']['control']:.2f}%. This is well above 50%, indicating strong attraction.")
    print(f" - White noise: {averages['White noise']['control']:.2f}%. This is below 50%, suggesting aversion.")
    print("Conclusion: The strong attraction to the tropical estuarine soundscape suggests it is the natural habitat.\n")
    
    # Analysis of CO2 impact
    print("Analysis of CO2 Impact:")
    control_tropical_avg = averages['Tropical estuarine soundscape']['control']
    co2_tropical_avg = averages['Tropical estuarine soundscape']['Co2']
    print(f" - For the Tropical estuarine soundscape, the behavior flips from attraction ({control_tropical_avg:.2f}%) to aversion ({co2_tropical_avg:.2f}%) under high CO2.")
    print("Conclusion: Elevated CO2 severely disturbs the larvae's ability to navigate to their natural habitat, which would make settlement less efficient.\n")

    print("Step 3: Evaluating the Answer Choices.\n")

    print("A: Incorrect. The data does not strongly support temperate reef as a natural habitat, and settlement behavior in the temperate reef is altered, not 'as efficiently'.")
    print("B: Incorrect. It correctly states settlement in the temperate reef is disturbed, but wrongly identifies the temperate reef as the primary natural habitat.")
    print("C: This option states that a natural habitat is the tropical estuarine and that at high CO2, the fish will not settle there efficiently. This aligns with our key findings. While calling the temperate reef a 'natural habitat' is weak, this option correctly identifies the tropical estuarine and the negative impact of CO2 on settlement there. This appears to be the best fit.")
    print("D: Correct, but too general. It doesn't capture the specific findings about habitat selection which is the focus of the study.")
    print("E: Incorrect. CO2 *does* disturb settlement in the tropical estuarine, and the habitat is likely not the temperate reef.")
    print("F: Incorrect. It correctly states the disturbance in the tropical estuarine but wrongly identifies the temperate reef as the natural habitat.")
    print("G: Incorrect, as C is a plausible answer.")
    print("H: Incorrect. The conclusion is correct but general, and the premise about habitats is weak and the same as C's. C offers a more specific and relevant conclusion about the primary habitat.")

    print("\nFinal Decision: Option C is the most accurate choice. It correctly identifies the tropical estuarine soundscape as a habitat and accurately describes the main finding of the experiment: that high CO2 levels will prevent the larvae from settling there efficiently.")
    
    print("<<<C>>>")

solve_fish_behavior_problem()