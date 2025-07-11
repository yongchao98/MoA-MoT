import math

# Based on the analytical derivation, the value of the integral is:
# I = 5^(1/4) / 12

# Let's define the components of the expression.
base = 12
power = 4
integral_value_numerator_base = 5

# Calculate the value of the integral I
integral_value = (integral_value_numerator_base**(1/power)) / base

# The full expression is (base^power) * (integral_value^power)
final_result = (base**power) * (integral_value**power)

# We print the final equation with each number.
# The float representation of 5^(1/4) is used for display purposes.
# The calculation of `final_result` is precise due to cancellation.

print(f"The equation to compute is: ({base})^{power} * (integral)^power")
print(f"After simplification, the integral evaluates to: ({integral_value_numerator_base}^(1/{power})) / {base}")
print("The final computation is:")
print(f"({base})^{power} * (({integral_value_numerator_base**(1/power)}) / {base})^{power} = {final_result:.0f}")
