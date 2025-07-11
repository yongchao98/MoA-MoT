# Define the constants from the derived formula for y(0)
# The formula is y(0) = c1 * (c2 ** exponent)
c1 = 3/4
c2 = 3/2
exponent = 6/5

# Calculate the final value of y(0)
y_at_0 = c1 * (c2 ** exponent)

# Print the final equation with its components and the result
print(f"The membrane's deflection at x=0, which is y(0), is calculated using the formula:")
print(f"y(0) = c1 * (c2 ** exponent)")
print(f"where c1 = {c1}, c2 = {c2}, and the exponent is {exponent}.")
print(f"Plugging in the numbers:")
print(f"y(0) = {c1} * ({c2} ** {exponent})")
print(f"y(0) = {y_at_0}")