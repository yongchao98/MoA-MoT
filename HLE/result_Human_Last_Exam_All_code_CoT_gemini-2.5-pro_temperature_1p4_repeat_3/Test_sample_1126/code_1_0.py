import numpy as np

def analyze_fish_data():
    """
    Analyzes fish larvae data to determine the effect of CO2 on habitat selection.
    """
    # Store the data in a dictionary
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

    print("Step 1: Identifying the Natural Habitat using Control data (current CO2 level)")
    print("-" * 70)
    
    # Calculate and print average for Tropical Estuarine (Control)
    tropical_control_data = data["Tropical estuarine soundscape"]["control"]
    tropical_control_avg = np.mean(tropical_control_data)
    # The 'f-string' formats the output to join the numbers with '+'
    print(f"Average time near Tropical Estuarine (Control): ({' + '.join(map(str, tropical_control_data))}) / {len(tropical_control_data)} = {tropical_control_avg:.2f}%")
    
    # Calculate and print average for Temperate Reef (Control)
    temperate_control_data = data["Temperate reef soundscape"]["control"]
    temperate_control_avg = np.mean(temperate_control_data)
    print(f"Average time near Temperate Reef (Control): ({' + '.join(map(str, temperate_control_data))}) / {len(temperate_control_data)} = {temperate_control_avg:.2f}%")
    
    print("\nConclusion 1: Under control conditions, larvae show a strong preference for the 'Tropical estuarine soundscape' (average > 50%), indicating it is their primary natural habitat. The response to the 'Temperate reef soundscape' is neutral (average is close to 50%).")

    print("\nStep 2: Assessing the impact of elevated CO2 on settlement behavior")
    print("-" * 70)
    
    # Calculate and print average for Tropical Estuarine (CO2)
    tropical_co2_data = data["Tropical estuarine soundscape"]["co2"]
    tropical_co2_avg = np.mean(tropical_co2_data)
    print(f"Average time near Tropical Estuarine (CO2): ({' + '.join(map(str, tropical_co2_data))}) / {len(tropical_co2_data)} = {tropical_co2_avg:.2f}%")
    
    print(f"\nConclusion 2: For the Tropical Estuarine habitat, the preference dropped from {tropical_control_avg:.2f}% to {tropical_co2_avg:.2f}%. This significant drop shows that elevated CO2 levels will impair the larvae's ability to settle in their natural habitat.")

    print("\nStep 3: Evaluating Answer Choices")
    print("-" * 70)
    print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
    print("\nAnalysis of Choice C:")
    print(" - The first part of the statement ('natural habitats are... tropical estuarine and the temperate reef') is an acceptable premise for the question, though the data shows a much stronger preference for the tropical estuarine.")
    print(f" - The second part ('will not settle in the tropical estuarine efficiently') is strongly supported by our analysis. The preference dropped from {tropical_control_avg:.2f}% to {tropical_co2_avg:.2f}%, which is a clear impairment of settlement efficiency.")
    print(" - Compared to other options, which contain clear contradictions with the data (e.g., claiming no disturbance or misidentifying the primary habitat), Choice C provides the most accurate conclusion based on the experimental results.")

analyze_fish_data()
<<<C>>>