def analyze_fish_larvae_data():
    """
    Calculates and prints the analysis of fish larvae behavior data
    under control and elevated CO2 conditions.
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
        print(f"Analysis for: {soundscape}")
        
        # Control condition
        control_data = conditions["control"]
        control_sum = sum(control_data)
        control_len = len(control_data)
        control_avg = control_sum / control_len
        control_eq = f"({ ' + '.join(map(str, control_data)) }) / {control_len}"
        print(f"  Control Average %: {control_eq} = {control_avg:.2f}%")

        # CO2 condition
        co2_data = conditions["co2"]
        co2_sum = sum(co2_data)
        co2_len = len(co2_data)
        co2_avg = co2_sum / co2_len
        co2_eq = f"({ ' + '.join(map(str, co2_data)) }) / {co2_len}"
        print(f"  CO2 Average %: {co2_eq} = {co2_avg:.2f}%")
        print("-" * 20)

    print("\n--- Interpretation ---")
    print("1. Natural Habitat Identification (based on Control conditions):")
    print("   - Tropical Estuarine (56.33%): Larvae show clear attraction. This is a likely natural habitat.")
    print("   - Temperate Reef (48.17%): Larvae show a neutral response (no attraction). This is likely not a preferred habitat.")
    print("   - White Noise (46.17%): Larvae show slight avoidance/neutrality as expected for an irrelevant sound.")
    
    print("\n2. Impact of Elevated CO2:")
    print("   - For the Tropical Estuarine soundscape, the larvae's behavior flips from strong attraction (56.33%) to avoidance (43.67%). This indicates that elevated CO2 significantly disturbs their ability to settle in their natural habitat.")
    print("   - For the Temperate Reef soundscape, the response remains neutral (48.17% vs 51.83%).")
    print("   - For White Noise, the behavior changes from neutral/avoidance to attraction (46.17% vs 53.83%), suggesting CO2 causes confusion.")

    print("\n--- Conclusion ---")
    print("The data shows that the natural habitat is the tropical estuarine. Elevated CO2 levels impair the larvae's ability to navigate to this habitat, making settlement less efficient. Option C best captures this critical finding, even if its premise about both locations being natural habitats is questionable. It correctly identifies that settlement in the tropical estuarine will be less efficient under future CO2 levels.")

    final_answer = "C"
    print(f"\nFinal Answer: <<<C>>>")

analyze_fish_larvae_data()