# Define the number of points of each color for a maximal configuration.
# This is based on the geometric construction that allows n=8.
num_red_points = 4
num_green_points = 2
num_yellow_points = 2

# Calculate the maximum value of n
max_n = num_red_points + num_green_points + num_yellow_points

# The problem asks to output the equation as well.
print(f"A valid configuration is r={num_red_points}, g={num_green_points}, y={num_yellow_points}.")
print(f"The maximum value of n is {num_red_points} + {num_green_points} + {num_yellow_points} = {max_n}.")