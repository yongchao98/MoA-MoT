# Define the numbers in the final expression for y(0)
a, b = 3, 4
c, d = 9, 4
e, f = 3, 5

# Calculate the value of y(0)
# y(0) = (3/4) * (9/4)^(3/5)
y0 = (a / b) * (c / d)**(e / f)

# Print the final equation with the numbers
print(f"The equation for the deflection at x=0 is:")
print(f"y(0) = ({a}/{b}) * ({c}/{d})**({e}/{f})")

# Print the calculated result
print("\nThe calculated value for y(0) is:")
print(y0)