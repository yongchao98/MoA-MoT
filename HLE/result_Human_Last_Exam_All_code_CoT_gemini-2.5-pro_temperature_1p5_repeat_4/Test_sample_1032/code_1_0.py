import math

# Step 1: Define the rope length and house area
rope_length = 7 / 2
house_area = 3.0

# Step 2: Calculate the area of the diamond shape reachable without obstacles
# The area of the region |x| + |y| <= L is 2 * L^2
diamond_area = 2 * (rope_length ** 2)

# Step 3: Explain the logic
# The key insight is that for taxi-cab distance, the length of the shortest path from the origin
# around a grid-aligned obstacle is the same as the direct taxi-cab distance |x| + |y|.
# This means a point is reachable if it's within the rope's length and not inside the house.
# The entire house is within the diamond since its furthest corner from the origin is (-2,-1) or (-1,-2),
# where the taxi-cab distance is |-2|+|-1| = 3, which is less than 3.5.

# Step 4: Calculate the final area
# The reachable area is the area of the diamond minus the area of the house.
final_area = diamond_area - house_area

# Print the explanation and the final calculation
print("The total area the horse could reach without any obstacles is a diamond shape.")
print(f"The taxi-cab rope length is L = {rope_length}.")
print(f"The area of this diamond is 2 * L^2 = 2 * {rope_length}^2 = {diamond_area}.")
print(f"The area of the house (obstacle) is {house_area} square units.")
print("Because the path length around the house corners is the same as the direct taxi-cab distance,")
print("the reachable area is simply the diamond's area minus the house's area.")
print("\nFinal Calculation:")
print(f"Reachable Area = {diamond_area} - {house_area} = {final_area}")
