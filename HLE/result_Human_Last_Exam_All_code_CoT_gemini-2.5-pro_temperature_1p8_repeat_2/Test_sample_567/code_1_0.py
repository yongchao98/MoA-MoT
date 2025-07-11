import math

# This script calculates the value 'a' where the volume constraint becomes the
# only obstruction for symplectically embedding the ellipsoid E(1,a) into a ball.
# This value is known to be the fourth power of the golden ratio, tau.

print("The problem asks for the value 'a' where the embedding capacity c(a) is determined purely by the volume constraint.")
print("This occurs when c(a) = sqrt(a).")
print("From the work of McDuff and Schlenk, this transition point is a = tau^4, where tau is the golden ratio.")
print("")

# The value can be expressed exactly as: a = (7 + 3*sqrt(5)) / 2
# We will define the numbers in this final equation.
seven = 7
three = 3
five = 5
two = 2

# Calculate the value of sqrt(5) for the equation
sqrt_five_val = math.sqrt(five)

# Calculate the numerator: 7 + 3*sqrt(5)
numerator_val = seven + three * sqrt_five_val

# Calculate the final value of 'a'
a_val = numerator_val / two

# Print the final equation with each of its constituent numbers displayed, as requested.
print(f"The final equation for 'a' is:")
print(f"a = ({seven} + {three} * sqrt({five})) / {two}")
print("")

# Print the step-by-step numerical calculation.
print(f"To find the numerical value, we calculate:")
print(f"a = ({seven} + {three} * {sqrt_five_val:.6f}) / {two}")
print(f"a = ({numerator_val:.6f}) / {two}")
print(f"a = {a_val:.6f}")
print("")

print(f"Thus, the value at which the only obstruction becomes the volume constraint is a = {a_val:.6f}")
