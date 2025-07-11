import math

# Earth's axial tilt in degrees, given in the problem.
epsilon_deg = 23.5

# As derived from the problem's conditions, the angular distance (d_sigma)
# between the two stars can be found using the formula: d_sigma = 180 - 2 * epsilon.
# Here we perform this calculation.

# The value of 2 * epsilon
two_epsilon = 2 * epsilon_deg

# The final angular distance
angular_distance = 180 - two_epsilon

# The problem asks to output each number in the final equation.
# The following print statement constructs and displays this equation along with the result.
print("The final equation for the angular distance is:")
# Using an f-string to embed the variables into the printed text.
# The numbers in the equation are 180, 2, the value of epsilon, and the final result.
print(f"180 - 2 * {epsilon_deg} = {angular_distance}")