import math

# The analytical solution leads to a parametric form for y(x).
# We found that to find y(0), we can use the formula:
# y(0) = (3/4) * (9/4)^(3/5)
# Here are the numbers in the final equation:
num1 = 3
den1 = 4
base_num = 9
base_den = 4
exp_num = 3
exp_den = 5

# Calculate the value
y0 = (num1 / den1) * (base_num / base_den)**(exp_num / exp_den)

# Print the final equation and the result
print("The deflection at x = 0 is given by the equation:")
# The instruction requested to output each number in the final equation.
# Here we display the equation structure with its components.
print(f"y(0) = ({num1}/{den1}) * ({base_num}/{base_den})^({exp_num}/{exp_den})")
print(f"y(0) = {y0}")
