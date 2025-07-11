import sys
import io

# Helper function to print to a string
def print_to_string(*args, **kwargs):
    old_stdout = sys.stdout
    sys.stdout = new_stdout = io.StringIO()
    print(*args, **kwargs)
    sys.stdout = old_stdout
    return new_stdout.getvalue()

# Data from the experiment
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

output = []
output.append("Step 1: Calculating Average Attraction for Each Scenario\n")

# Calculate and store averages
averages = {}
for soundscape, conditions in data.items():
    control_avg = sum(conditions['control']) / len(conditions['control'])
    co2_avg = sum(conditions['co2']) / len(conditions['co2'])
    averages[soundscape] = {"control": control_avg, "co2": co2_avg}
    
    output.append(f"--- {soundscape} ---")
    control_data_str = ', '.join(map(str, conditions['control']))
    co2_data_str = ', '.join(map(str, conditions['co2']))
    
    # Generate the equation string for control average
    control_equation_str = f"({ ' + '.join(map(str, conditions['control'])) }) / {len(conditions['control'])}"
    output.append(print_to_string(f"Control Average = {control_equation_str} = {control_avg:.2f}%"))
    
    # Generate the equation string for CO2 average
    co2_equation_str = f"({ ' + '.join(map(str, conditions['co2'])) }) / {len(conditions['co2'])}"
    output.append(print_to_string(f"CO2 Average     = {co2_equation_str} = {co2_avg:.2f}%"))


output.append("\nStep 2: Analysis and Interpretation\n")
# Interpretation of natural habitat
tropical_control_avg = averages["Tropical estuarine soundscape"]["control"]
temperate_control_avg = averages["Temperate reef soundscape"]["control"]

output.append(f"1. Natural Habitat Identification:")
output.append(f"   - Under control conditions, larvae showed strong attraction to the 'Tropical estuarine soundscape' (average: {tropical_control_avg:.2f}%). This indicates it is a natural habitat.")
output.append(f"   - For the 'Temperate reef soundscape', the larvae's response was neutral (average: {temperate_control_avg:.2f}%, close to 50%). This indicates it is likely NOT a preferred natural habitat.")

# Interpretation of CO2 effect
tropical_co2_avg = averages["Tropical estuarine soundscape"]["co2"]
output.append(f"2. Impact of Elevated CO2:")
output.append(f"   - In the 'Tropical estuarine soundscape', elevated CO2 caused a dramatic shift from attraction ({tropical_control_avg:.2f}%) to aversion ({tropical_co2_avg:.2f}%).")
output.append(f"   - This means the larvae will NOT settle efficiently in their natural habitat under future CO2 levels.")

output.append("\nStep 3: Evaluating Answer Choices\n")
output.append("Based on the analysis:")
output.append(" - Choice A is incorrect because the temperate reef is not indicated as a natural habitat.")
output.append(" - Choice B is incorrect because the temperate reef is not indicated as a natural habitat.")
output.append(" - Choice C correctly states that fish will not settle efficiently in the tropical estuarine under high CO2. Although its premise that the temperate reef is also a habitat is not supported by the data, this statement captures the most significant conclusion from the experiment.")
output.append(" - Choice D is true but too general.")
output.append(" - Choices E and F are incorrect because they state the temperate reef is the natural habitat.")
output.append(" - Choice H is incorrect because the temperate reef is not indicated as a natural habitat.")
output.append("\nConclusion: Choice C provides the most accurate and specific conclusion about the key findings of the experiment, despite a questionable premise about one of the habitats.")

print(''.join(output))

print("\n<<<C>>>")