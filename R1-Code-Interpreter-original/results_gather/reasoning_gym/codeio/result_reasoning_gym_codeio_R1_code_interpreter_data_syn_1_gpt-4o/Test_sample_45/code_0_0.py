import numpy as np
from math import atan, sin, sqrt

def gkk(theta, units='radians'):
    # If theta is a list, makes it into a numpy array
    if isinstance(theta, list):
        theta = np.array(theta)
    
    # Calculates a constant for easy reference
    at2 = atan(1.0 / sqrt(2.0))
    
    # Splits calculation based on how theta is provided; does here if it's a numpy array
    gkk_fracq = np.piecewise(theta, [np.abs(theta) <= at2, np.abs(theta) > at2], \
                            [-1.0, lambda x: (3.0 * (np.sin(x))**2 - 2.0) / (3.0 * (np.sin(x))**2)])
    gkk_evpa = 0.5 * np.arctan2(0.0, gkk_fracq)
    
    # Return fracq and evpa
    return gkk_fracq, gkk_evpa

# main function
def main_solution(theta, units):
    # Convert theta to a list if it's a single value
    if not isinstance(theta, list):
        theta = [theta]
    
    # Call the gkk function
    gkk_fracq, gkk_evpa = gkk(theta, units)
    
    # Convert the results to JSON serializable format
    gkk_fracq = gkk_fracq.tolist() if isinstance(gkk_fracq, np.ndarray) else gkk_fracq
    gkk_evpa = gkk_evpa.tolist() if isinstance(gkk_evpa, np.ndarray) else gkk_evpa
    
    # Return the results as a dictionary
    return {"gkk_fracq": gkk_fracq, "gkk_evpa": gkk_evpa}

# Given input
input_data = {'theta': [1.5110653204394262, 5.168270952748489], 'units': 'radians'}

# Calculate the output
output = main_solution(input_data['theta'], input_data['units'])

# Print the output
print(output)