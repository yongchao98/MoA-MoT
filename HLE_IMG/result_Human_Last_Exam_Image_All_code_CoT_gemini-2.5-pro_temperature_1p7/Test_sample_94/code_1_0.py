import math

# Step 1: Analyze the Monte Carlo simulation data
# AGI reported the following number of points falling into circles out of 1000 trials:
points_in_circles = [730, 740, 735, 732, 739]
total_points = 1000
# Calculate the average number of points
avg_points_in_circles = sum(points_in_circles) / len(points_in_circles)
# This gives the ratio of the area of the circles to the total area
ratio_circles_total = avg_points_in_circles / total_points
print(f"Step 1: The average number of points in circles is ({points_in_circles[0]} + {points_in_circles[1]} + {points_in_circles[2]} + {points_in_circles[3]} + {points_in_circles[4]}) / 5 = {avg_points_in_circles}.")
print(f"This gives the ratio of the circle area to total area: {avg_points_in_circles} / {total_points} = {ratio_circles_total:.4f}.")
print("-" * 30)

# Step 2: Use the information about the green rectangles
# The area of one green rectangle is ~4% of the total area.
ratio_green_total = 0.04
print(f"Step 2: The ratio of a green rectangle's area to the total area is given as {ratio_green_total}.")
print("-" * 30)

# Step 3: Establish a relationship between the component dimensions
# From the image, there are 9 circles. Total circle area = 9 * pi * r^2
# From the image, h_g = 2r. Area of one green rectangle = w_g * 2r.
# The ratio of these areas is independent of the unknown total area:
# Area_green / Area_circles = (w_g * 2r) / (9 * pi * r^2) = 2 * w_g / (9 * pi * r)
# This must equal the ratio of their percentages: ratio_green_total / ratio_circles_total
ratio_green_circles = ratio_green_total / ratio_circles_total
# So, we can solve for w_g / r:
# w_g / r = (ratio_green_circles) * (9 * pi / 2)
w_div_r_ratio = ratio_green_circles * (9 * math.pi / 2)
print(f"Step 3: From the data, we derive a relationship between the rectangle width w_g and circle radius r:")
print(f"w_g / r = ({ratio_green_total} / {ratio_circles_total:.4f}) * (9 * pi / 2)")
print(f"w_g / r ≈ {w_div_r_ratio:.4f}")
print("-" * 30)

# Step 4: Find the integer values for r and w_g
# Since r and w_g are integers, we need to find a fraction that approximates 0.7692.
# The fraction 10/13 ≈ 0.76923 is an excellent match.
r = 13
w_g = 10
h_g = 2 * r
print(f"Step 4: By finding the best integer fraction for {w_div_r_ratio:.4f}, we find:")
print(f"Radius of circles, r = {r} cm")
print(f"Width of green rectangle, w_g = {w_g} cm")
print(f"Height of green rectangle, h_g = 2 * r = {h_g} cm")
print("-" * 30)

# Step 5: Calculate the total area of the outer rectangle
# Area of one green rectangle = w_g * h_g
area_green = w_g * h_g
# Total area = Area of green rectangle / ratio_green_total
total_area = area_green / ratio_green_total
print(f"Step 5: The area of one green rectangle is {w_g} * {h_g} = {area_green} sq cm.")
print(f"The total area of the outer rectangle is therefore {area_green} / {ratio_green_total} = {total_area} sq cm.")
print("-" * 30)

# Step 6: Determine the dimensions x and y of the outer rectangle
# From the layout, the width x is the span of 3 packed circles (6r) plus the width of one green rectangle (w_g).
x = 6 * r + w_g
# The height y can be calculated from the total area and the width x.
y = total_area / x
print("Step 6: Based on the visual arrangement:")
print(f"The total width x = 6 * r + w_g = 6 * {r} + {w_g} = {x} cm.")
print(f"The total height y = Total Area / x = {total_area} / {x} ≈ {y:.2f} cm.")
print("-" * 30)

# Step 7: Final Answer
# Round the dimensions to the nearest cm.
x_rounded = round(x)
y_rounded = round(y)
print("Step 7: The final answer, rounded to the nearest cm, is...")
final_answer = f"{x_rounded}:{y_rounded}"
print(final_answer)