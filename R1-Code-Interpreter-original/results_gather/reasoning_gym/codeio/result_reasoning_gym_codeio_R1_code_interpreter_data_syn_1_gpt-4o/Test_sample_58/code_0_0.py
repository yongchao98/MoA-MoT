import numpy as np

def convert_measurement(pivot, end_position, value, isRotational, toRotational):
    # Convert lists to numpy arrays for vector operations
    pivot = np.array(pivot)
    end_position = np.array(end_position)
    
    # Calculate the radius vector and its magnitude
    radius_vector = end_position - pivot
    radius_magnitude = np.linalg.norm(radius_vector)
    
    # Determine the converted value based on the desired mode
    if isRotational and not toRotational:
        # Rotational to Translational
        converted_value = value * radius_magnitude
    elif not isRotational and toRotational:
        # Translational to Rotational
        converted_value = value / radius_magnitude
    else:
        # No conversion needed
        converted_value = value
    
    # Return the result as a dictionary
    return {
        "converted_value": converted_value,
        "isRotational": toRotational
    }

# Given input
input_data = {
    'pivot': [42.7194000671706, -3.7942742415070683, 14.920994022681427],
    'end_position': [-38.03998584388315, 77.28626929945787, 60.99518157730225],
    'value': 57.67340626379919,
    'isRotational': False,
    'toRotational': False
}

# Calculate the output
output = convert_measurement(**input_data)
print(output)