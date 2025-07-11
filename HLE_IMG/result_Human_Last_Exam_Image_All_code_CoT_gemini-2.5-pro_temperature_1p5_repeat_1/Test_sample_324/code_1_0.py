import math

# Step 1: Assign the determined values to variables.
# n0 is the index of the base plot.
n0 = 1
# n1 is the index of the plot with reduced annulus pressure gradient.
n1 = 2
# n6 is the index of the plot with reduced tube-to-annulus viscosity ratio.
n6 = 6
# The ratio of maximum velocities under the specified conditions simplifies to ln(2).
v_ratio = math.log(2)

# Step 2: Calculate the components of the expression.
prefactor = (n1 / n6) ** n0
exp_term = math.exp(v_ratio)

# Step 3: Calculate the final result.
final_result = prefactor * exp_term

# Step 4: Print the final equation with each number and the result.
# The problem asks to output each number in the final equation.
# Here we print the equation in its symbolic and numerical forms.
print(f"Final Equation: (n1 / n6)^n0 * exp(v_a_max / v_t_max)")
print(f"Substituting values: ({n1} / {n6})^{n0} * exp(ln(2))")
print(f"Calculation: (1/3) * 2")
print(f"Result: {final_result}")
<<<0.6666666666666667>>>