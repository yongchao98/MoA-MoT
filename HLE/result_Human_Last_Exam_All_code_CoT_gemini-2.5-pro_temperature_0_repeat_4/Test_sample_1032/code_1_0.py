# Define the parameters of the problem
rope_length = 7 / 2
house_area = 3

# The area the horse can reach without any obstacles is a diamond shape.
# The area of a diamond in taxi-cab geometry is 2 * L^2, where L is the radius (rope length).
diamond_area = 2 * (rope_length**2)

# The problem states the horse cannot enter the house.
# Our analysis shows that because the rope is tied to a corner of the house,
# the shortest path to any point (x,y) is still |x|+|y|.
# Therefore, the reachable area is the area of the diamond minus the area of the house.
reachable_area = diamond_area - house_area

# Print the final equation with each number
print(f"The total reachable area is calculated by the area of the diamond minus the area of the house.")
print(f"Area_diamond = 2 * ({rope_length})**2 = {diamond_area}")
print(f"Area_house = {house_area}")
print(f"Final Area = {diamond_area} - {house_area} = {reachable_area}")
