import numpy as np

def analyze_settlement_data():
    """
    Analyzes fish larvae settlement data to determine the impact of CO2.
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
        control_data = conditions["control"]
        co2_data = conditions["co2"]

        # Calculate averages
        control_avg = np.mean(control_data)
        co2_avg = np.mean(co2_data)

        # Print results
        print(f"  Control Group Average: {control_avg:.2f}%")
        print(f"  High CO2 Group Average: {co2_avg:.2f}%")
        print("-" * 20)
    
    print("\n--- Conclusions from Data ---")
    print("1. Natural Habitat Identification (Control Conditions):")
    print("   - Tropical Estuarine: Average is 56.33%, with initial strong attraction (68%, 63%). This indicates it is a preferred natural habitat.")
    print("   - Temperate Reef: Average is 48.17%, which is close to 50% (indifference). Not a strongly preferred habitat based on this sound cue.")
    print("   - White Noise: Average is 46.17%, with initial strong avoidance (30%). This is an irrelevant/undesirable sound, as expected.")

    print("\n2. Impact of High CO2:")
    print("   - In the Tropical Estuarine (preferred) habitat, high CO2 caused a major shift from attraction (56.33%) to avoidance (43.67%). This is a clear disruption, making settlement inefficient.")
    print("   - For White Noise, high CO2 reversed the behavior from avoidance (46.17%) to attraction (53.83%). This confirms sensory function is impaired.")
    print("   - In the Temperate Reef, the change is from indifference (48.17%) to slight attraction (51.83%). Behavior is altered.")
    
    print("\n--- Evaluating the Best Answer ---")
    print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
    print("This aligns perfectly with our findings. The most significant result is the disruption of settlement in the clearly preferred tropical estuarine habitat.")

analyze_settlement_data()
<<<C>>>