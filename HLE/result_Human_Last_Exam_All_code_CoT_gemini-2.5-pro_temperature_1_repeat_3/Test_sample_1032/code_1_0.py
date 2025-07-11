# Define the length of the rope
rope_length = 7 / 2

# The total area the horse could reach without any obstacles is a diamond shape
# defined by |x| + |y| <= rope_length. The area of this diamond is 2 * L^2.
diamond_area = 2 * rope_length**2

# The house consists of three 1x1 unit squares, so its total area is 3.
house_area = 3

# As shown in the derivation, the shortest taxi-cab path length from the origin
# to any point (x,y) around the house is still |x| + |y|.
# Therefore, the reachable region is the area of the diamond minus the area of the house.
reachable_area = diamond_area - house_area

# Print the equation with all the numbers
print(f"The total reachable area is calculated by the equation:")
print(f"2 * ({rope_length})**2 - {house_area} = {reachable_area}")

# Print the final numerical answer
# print(f"The area of the region the horse can reach is: {reachable_area}")