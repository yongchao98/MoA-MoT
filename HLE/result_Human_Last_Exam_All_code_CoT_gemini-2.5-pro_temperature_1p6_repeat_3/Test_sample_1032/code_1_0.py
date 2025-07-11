# Define the taxi-cab length of the rope
rope_length = 7 / 2

# Define the area of the house obstacle
# The problem states it consists of three unit squares.
house_area = 3.0

# Calculate the total area reachable by the horse without any obstacles.
# This area is a diamond shape defined by |x| + |y| <= rope_length.
# The area of such a diamond is 2 * L^2.
diamond_area = 2 * rope_length**2

# The crucial insight is that for this specific geometry, the taxi-cab distance
# from the origin to any point (x, y) around the house obstacle is still |x| + |y|.
# Therefore, the horse can reach any point within the diamond |x| + |y| <= 3.5,
# as long as that point is not inside the house.
# The entire house is within this diamond, so we subtract its area.
final_area = diamond_area - house_area

# Output the steps of the calculation
print(f"The taxi-cab length of the rope is L = {rope_length}.")
print(f"The total area the horse could reach with no obstacles is a diamond with area 2 * L^2.")
print(f"Area of the diamond = 2 * ({rope_length})^2 = {diamond_area}.")
print(f"The area of the house (obstacle) is {house_area}.")
print(f"The final reachable area is the area of the diamond minus the area of the house.")
print(f"Final Area = {diamond_area} - {house_area} = {final_area}")
