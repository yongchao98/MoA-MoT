import numpy as np
from scipy.optimize import minimize

# Define the objective function to minimize
def objective(vars):
    x = vars[:4]
    y = vars[4:]
    norm_x = np.linalg.norm(x)
    norm_y = np.linalg.norm(y)
    norm_x_minus_y = np.linalg.norm(x - y)
    
    # Calculate the squared differences from the target norms
    diff_x = (norm_x - 7.055491021319009) ** 2
    diff_y = (norm_y - 8.225669343609868) ** 2
    diff_x_minus_y = (norm_x_minus_y - 13.830100438670305) ** 2
    
    return diff_x + diff_y + diff_x_minus_y

# Initial guess for the vectors
initial_guess = np.random.rand(8)

# Perform the optimization
result = minimize(objective, initial_guess, method='BFGS')

# Extract the vectors from the result
x_values = result.x[:4]
y_values = result.x[4:]

print(x_values, y_values)