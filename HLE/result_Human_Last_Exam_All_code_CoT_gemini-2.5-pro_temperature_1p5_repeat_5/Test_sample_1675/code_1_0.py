# Define the sets of colored points for a valid n=8 configuration.
# This configuration has |R|=3, |G|=4, |Y|=1.

R = {'r1': (0, 10), 'r2': (10, -5), 'r3': (-10, -5)}
G = {'g1': (0, 1), 'g2': (3, 0), 'g3': (0, -1), 'g4': (-3, 0)}
Y = {'y1': (0, 0)}

# Get the number of points of each color
num_R = len(R)
num_G = len(G)
num_Y = len(Y)

# Calculate the total number of points n
n = num_R + num_G + num_Y

# Print the breakdown of the point distribution
print("A valid configuration exists with:")
print(f"Number of red points: {num_R}")
print(f"Number of green points: {num_G}")
print(f"Number of yellow points: {num_Y}")
print("")

# Print the equation for the total number of points
print("The maximum value of n is the sum of these numbers:")
print(f"{num_R} + {num_G} + {num_Y} = {n}")