import numpy as a

def calculate_and_print_average(name, data_dict):
    """Calculates and prints the average for control and CO2 conditions."""
    print(f"--- Analysis for {name} ---")
    
    control_data = data_dict['control']
    co2_data = data_dict['co2']
    
    # Calculate averages
    avg_control = a.mean(control_data)
    avg_co2 = a.mean(co2_data)
    
    # Format numbers for the equation string
    control_nums_str = ' + '.join(map(str, control_data))
    co2_nums_str = ' + '.join(map(str, co2_data))
    
    # Print the detailed calculation for control
    print("Control Condition:")
    print(f"Equation: ({control_nums_str}) / {len(control_data)}")
    print(f"Result: Average time spent near speaker = {avg_control:.2f}%")
    if avg_control > 50:
        print("Conclusion: Larvae show attraction, suggesting this is a suitable habitat.\n")
    else:
        print("Conclusion: Larvae do not show attraction.\n")

    # Print the detailed calculation for CO2
    print("High CO2 Condition:")
    print(f"Equation: ({co2_nums_str}) / {len(co2_data)}")
    print(f"Result: Average time spent near speaker = {avg_co2:.2f}%")
    if avg_control > 50 and avg_co2 < avg_control:
        print("Conclusion: Attraction is reduced or lost, suggesting settlement will be less efficient.\n")
    else:
        print("Conclusion: Behavior is altered by high CO2.\n")


# Data from the problem description
# Note: Corrected a typo from 'Day2' to 'Day21' for Temperate Reef
data = {
    "Temperate Reef Soundscape": {
        "control": [40, 50, 49, 53, 49, 48],
        "co2": [60, 50, 51, 47, 51, 52]
    },
    "Tropical Estuarine Soundscape": {
        "control": [68, 63, 55, 51, 49, 52],
        "co2": [32, 37, 45, 49, 51, 48]
    }
}

# Perform and print the analysis
calculate_and_print_average("Temperate Reef Soundscape", data["Temperate Reef Soundscape"])
calculate_and_print_average("Tropical Estuarine Soundscape", data["Tropical Estuarine Soundscape"])

print("--- Final Answer Rationale ---")
print("The calculations show that under control conditions, larvae are attracted to the Tropical Estuarine sound (avg > 50%), but not the Temperate Reef sound (avg < 50%).")
print("Under high CO2, the attraction to the Tropical Estuarine sound is lost (avg drops from 56.33% to 43.67%).")
print("This directly supports the conclusion in option C that at high CO2 levels, 'the fish will not settle in the tropical estuarine efficiently as in the year 2024'.")
print("While the premise in option C about the temperate reef also being a habitat is not supported by the data, the conclusion about the tropical estuarine is the most significant and accurate finding.")

<<<C>>>