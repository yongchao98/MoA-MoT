# The equation describing the deflection at x=0 is derived as y(0) = (3/4) * (9/4)^(3/5).
# Here we define the numbers in this final equation.
num1 = 3
den1 = 4
num2 = 9
den2 = 4
num_pow = 3
den_pow = 5

# We output the numbers that form the final equation as requested.
print(f"The derived equation for the deflection at x=0 is y(0) = ({num1}/{den1}) * ({num2}/{den2})^({num_pow}/{den_pow})")

# Now, we calculate the numerical value.
value = (num1 / den1) * (num2 / den2)**(num_pow / den_pow)

# Finally, we print the result for y(0).
print(f"\nThe membrane's deflection at x = 0 is y(0) = {value}")