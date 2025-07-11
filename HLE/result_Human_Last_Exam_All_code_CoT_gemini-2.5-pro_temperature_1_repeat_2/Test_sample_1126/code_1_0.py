def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the effect of CO2 on habitat selection.
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

    for soundscape, conditions in data.items():
        print(f"Analyzing: {soundscape}")
        
        # Control group analysis
        control_data = conditions["control"]
        control_sum = sum(control_data)
        control_len = len(control_data)
        control_avg = control_sum / control_len
        control_sum_str = " + ".join(map(str, control_data))
        print(f"Control Average: ({control_sum_str}) / {control_len} = {control_avg:.2f}%")

        # CO2 group analysis
        co2_data = conditions["co2"]
        co2_sum = sum(co2_data)
        co2_len = len(co2_data)
        co2_avg = co2_sum / co2_len
        co2_sum_str = " + ".join(map(str, co2_data))
        print(f"CO2 Average:     ({co2_sum_str}) / {co2_len} = {co2_avg:.2f}%\n")

    print("--- Conclusion ---")
    print("1. Natural Habitat: Under control conditions, larvae show significant attraction only to the 'Tropical estuarine soundscape' (average > 50%), indicating it is their natural habitat.")
    print("2. Effect of CO2: In the 'Tropical estuarine soundscape', high CO2 levels cause the average time spent near the speaker to drop from 56.33% to 43.67%. This demonstrates a loss of attraction.")
    print("3. Final Answer: This loss of attraction means the larvae will not settle efficiently in their natural habitat. Choice C accurately describes this outcome, stating that 'the fish will not settle in the tropical estuarine efficiently'.")

analyze_fish_behavior()
print("<<<C>>>")