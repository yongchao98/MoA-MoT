import math

# Step 1: Define the rope length and calculate the unconstrained area.
L = 7/2
unconstrained_area = 2 * L**2

# Step 2: Define the area of the house.
# The house is a 2x2 square with a 1x1 square removed from a corner.
# Area = 2*2 - 1*1 = 3
house_area = 3.0

# Step 3: Analyze the effect of the obstacle.
# As explained in the plan, the re-entrant corner at (-1,-1) allows paths to be
# rerouted without increasing the taxi-cab distance.
# So, the reachable area is the unconstrained area minus the house area.

# Step 4: Calculate the final area.
final_area = unconstrained_area - house_area

# Step 5: Print the final equation with all numbers.
print("The taxi-cab rope length is L = 7/2.")
print("The unconstrained reachable area is a diamond shape with area given by the formula 2 * L^2.")
print("The area of the house is 3.")
print("Because of the house's shape, the taxi-cab distance to any point outside the house is not increased.")
print("The final reachable area is the unconstrained area minus the house area.")
print(f"Final Area = 2 * (7/2)^2 - 3 = 2 * {L**2} - {house_area} = {unconstrained_area} - {house_area} = {final_area}")
