# Based on the derivation, the original expression simplifies to:
# (12)^4 * ( integral_value )^4
# where integral_value = (5^(1/4)) / 12.

# Let's define the numbers from the simplified equation.
# Expression is base^power * ((num_base^(1/4))/den)^power
base = 12
power = 4
num_base = 5
den = 12

# The calculation simplifies to base^power * (num_base / den^power)
# which is 12^4 * (5 / 12^4), which simplifies to 5.
result = base**power * (num_base / (den**power))

# We print the final equation with the numbers involved.
print(f"The calculation is: ({base})^{power} * (({num_base}^(1/4)) / {den})^{power}")
print(f"This simplifies to: ({base})^{power} * {num_base} / ({den})^{power} = {int(result)}")