import numpy as np

def analyze_fish_larvae_data():
    """
    Analyzes fish larvae settlement data under different CO2 levels and soundscapes.
    """

    # 1. Organize the data
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

    print("Step 1: Analyzing behavior under control (current) conditions to identify natural habitat.\n")
    print("A value > 50% suggests attraction, a value < 50% suggests avoidance.")

    # 2. Calculate and print control averages
    tropical_control_avg = np.mean(data["Tropical estuarine soundscape"]["control"])
    print(f"Tropical Estuarine (Control): Data = {data['Tropical estuarine soundscape']['control']}")
    print(f"Average time spent near Tropical Estuarine soundscape in control: {tropical_control_avg:.2f}%")
    print("Conclusion: Larvae show a clear attraction to the Tropical Estuarine soundscape.\n")

    temperate_control_avg = np.mean(data["Temperate reef soundscape"]["control"])
    print(f"Temperate Reef (Control): Data = {data['Temperate reef soundscape']['control']}")
    print(f"Average time spent near Temperate Reef soundscape in control: {temperate_control_avg:.2f}%")
    print("Conclusion: Larvae show no significant preference for the Temperate Reef soundscape (value is close to 50%).\n")
    
    print("--> Based on control data, the Tropical Estuarine soundscape is a natural habitat.\n")
    
    print("Step 2: Analyzing behavior under elevated CO2 conditions to see the effect.\n")

    # 3. Calculate and print CO2 averages and compare
    tropical_co2_avg = np.mean(data["Tropical estuarine soundscape"]["co2"])
    print(f"Tropical Estuarine (CO2): Data = {data['Tropical estuarine soundscape']['co2']}")
    print(f"Average time spent near Tropical Estuarine soundscape with elevated CO2: {tropical_co2_avg:.2f}%")
    print(f"Comparison: The attraction seen in control conditions (avg {tropical_control_avg:.2f}%) is lost and turns into avoidance (avg {tropical_co2_avg:.2f}%) under high CO2.\n")
    
    print("Step 3: Evaluating the Answer Choices.\n")
    print("A: False. Behavior in temperate reef is altered, not 'as efficiently'.")
    print("B: False. States temperate reef is the natural habitat, which contradicts the strong attraction shown for tropical estuarine.")
    print("C: True. This choice correctly states that fish will not settle efficiently in the tropical estuarine habitat under future CO2 levels. Our analysis shows attraction (efficient settlement) is lost and reversed.")
    print("D: True, but too general. C is more specific and accurately describes the key finding.")
    print("E: False. Settlement in tropical estuarine is clearly disturbed.")
    print("F: False. States temperate reef is the natural habitat.")
    print("H: False. The claim that both are habitats is not strongly supported. C is more precise about the demonstrated effect.")
    
    print("\nFinal conclusion is that elevated CO2 disrupts the larvae's ability to use sound cues for settlement, specifically by making them unable to identify their natural habitat (Tropical Estuarine).")
    
analyze_fish_larvae_data()
<<<C>>>