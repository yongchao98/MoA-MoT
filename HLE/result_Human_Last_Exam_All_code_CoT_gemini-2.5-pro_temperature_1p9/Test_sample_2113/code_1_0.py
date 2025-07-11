import numpy as np

def z1_func(t):
  """Calculates the trajectory z1(t) = cos(t)"""
  return np.cos(t)

def y_func(t):
  """Calculates the entanglement echo y(t) = sin(2t) / (cos(2t))^(3/2)"""
  return np.sin(2 * t) / (np.cos(2 * t))**1.5

# Time at which to evaluate the functions
t_val = np.pi / 8

# Calculate the values of z1 and y at t_val
z1_val = z1_func(t_val)
y_val = y_func(t_val)

# Calculate the final expression
result = (z1_val / y_val)**2

# Print the components of the final equation and the result
print(f"The final calculation is (z_1(π/8) / y(π/8))^2.")
print(f"Value of z_1(π/8) = cos(π/8) ≈ {z1_val:.4f}")
print(f"Value of y(π/8) = 2^(1/4) ≈ {y_val:.4f}")
print(f"Thus, we compute ({z1_val:.4f} / {y_val:.4f})^2 ≈ {result:.4f}")

# The exact analytical result is (sqrt(2)+1)/4
exact_result_str = "(√2 + 1) / 4"
print(f"\nThe exact analytical result is {exact_result_str}.")
print(f"The final numerical value is: {result}")