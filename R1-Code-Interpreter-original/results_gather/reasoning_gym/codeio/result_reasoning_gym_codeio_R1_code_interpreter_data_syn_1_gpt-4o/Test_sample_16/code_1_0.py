import numpy as np

# Vectors obtained from the optimization
x_values = [5.38827146, 2.59932879, 0.19706185, 3.73511845]
y_values = [-2.34057426, -3.41981678, 5.56316317, -4.42034236]

# Calculate norms
norm_x = np.linalg.norm(x_values)
norm_y = np.linalg.norm(y_values)
norm_x_minus_y = np.linalg.norm(np.array(x_values) - np.array(y_values))

# Print the calculated norms
print(norm_x, norm_y, norm_x_minus_y)