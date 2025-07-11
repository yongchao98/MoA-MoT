def analyze_fish_data():
    """
    Analyzes fish larvae behavior data to determine habitat preference and the effect of CO2.
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

    # Step 1: Analyze Control Conditions to find Natural Habitat
    print("Step 1: Analyzing behavior under 'Control' (current CO2) conditions.")
    
    # Tropical Estuarine (Control)
    tropical_control_vals = data["Tropical estuarine soundscape"]["control"]
    tropical_control_avg = sum(tropical_control_vals) / len(tropical_control_vals)
    tropical_control_eq = f"({ ' + '.join(map(str, tropical_control_vals)) }) / {len(tropical_control_vals)}"
    print(f"Average time at Tropical Estuarine soundscape: {tropical_control_eq} = {tropical_control_avg:.2f}%")

    # Temperate Reef (Control)
    temperate_control_vals = data["Temperate reef soundscape"]["control"]
    temperate_control_avg = sum(temperate_control_vals) / len(temperate_control_vals)
    temperate_control_eq = f"({ ' + '.join(map(str, temperate_control_vals)) }) / {len(temperate_control_vals)}"
    print(f"Average time at Temperate Reef soundscape: {temperate_control_eq} = {temperate_control_avg:.2f}%")
    
    # White Noise (Control)
    noise_control_vals = data["White noise"]["control"]
    noise_control_avg = sum(noise_control_vals) / len(noise_control_vals)
    noise_control_eq = f"({ ' + '.join(map(str, noise_control_vals)) }) / {len(noise_control_vals)}"
    print(f"Average time at White Noise: {noise_control_eq} = {noise_control_avg:.2f}%")
    
    print("\nConclusion for Step 1: Larvae are strongly attracted to the Tropical Estuarine soundscape (>50%) and avoid/are indifferent to White Noise (<50%). This suggests the Tropical Estuarine is a preferred natural habitat.\n")
    
    # Step 2: Analyze the effect of elevated CO2
    print("---")
    print("Step 2: Analyzing behavior under elevated 'Co2' conditions.")

    # Tropical Estuarine (Co2)
    tropical_co2_vals = data["Tropical estuarine soundscape"]["Co2"]
    tropical_co2_avg = sum(tropical_co2_vals) / len(tropical_co2_vals)
    tropical_co2_eq = f"({ ' + '.join(map(str, tropical_co2_vals)) }) / {len(tropical_co2_vals)}"
    print(f"Average time at Tropical Estuarine soundscape: {tropical_co2_eq} = {tropical_co2_avg:.2f}%")
    print(f"Result: Attraction decreased from {tropical_control_avg:.2f}% to {tropical_co2_avg:.2f}%.\n")

    # White Noise (Co2)
    noise_co2_vals = data["White noise"]["Co2"]
    noise_co2_avg = sum(noise_co2_vals) / len(noise_co2_vals)
    noise_co2_eq = f"({ ' + '.join(map(str, noise_co2_vals)) }) / {len(noise_co2_vals)}"
    print(f"Average time at White Noise: {noise_co2_eq} = {noise_co2_avg:.2f}%")
    print(f"Result: Behavior changed from avoidance/indifference ({noise_control_avg:.2f}%) to attraction ({noise_co2_avg:.2f}%).\n")

    print("Conclusion for Step 2: Elevated CO2 impairs the larvae's sensory abilities. They are no longer attracted to their natural habitat sound and become attracted to irrelevant noise. This is a clear 'disturbance' of settlement behavior.\n")
    
    print("---")
    print("Final Conclusion:")
    print("The data shows that elevated CO2 will disturb the settlement of the examined fish. The disturbance is demonstrated by a loss of attraction to their likely natural habitat (Tropical Estuarine) and a new, abnormal attraction to white noise. Option H best summarizes these findings by stating that the settlement is disturbed and that both the tropical estuarine and temperate reef are potential habitats.")

# Execute the analysis
analyze_fish_data()