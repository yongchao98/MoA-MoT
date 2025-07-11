import math

# Define the given values
radius = 34
side_length = 17

# Step 1: Calculate theta in radians. tan(theta) = side_length / radius
theta = math.atan(side_length / radius)

# Step 2: Calculate the area of one of the shaded regions.
# Area of one shaded region = Area(Segment) + Area(Hat Triangle)
# Area(Segment) = Area(Sector) - Area(Inner Triangle)
# Area(Sector) has angle 2*theta
area_sector = radius**2 * theta

# Area(Inner Triangle OP0P2) = 1/2 * R^2 * sin(2*theta)
# Using sin(2*theta) = 2*sin(theta)*cos(theta). From tan(theta), we have sin(theta) and cos(theta).
# Alternatively, sin(2*theta) = 2*tan(theta)/(1+tan(theta)^2)
sin_2theta = (2 * (side_length / radius)) / (1 + (side_length / radius)**2)
area_inner_triangle = 0.5 * radius**2 * sin_2theta

# Area of Segment
area_segment = area_sector - area_inner_triangle

# Area(Hat Triangle P0T_bP2) = Area(Polygon OP0T_bP2) - Area(Inner Triangle OP0P2)
# Area(Polygon OP0T_bP2) is sum of Area(Triangle OP0T_b) and Area(Triangle OT_bP2)
# Both of these triangles have an area of 0.5 * radius * side_length = 289
area_polygon = 2 * (0.5 * radius * side_length)
area_hat_triangle = area_polygon - area_inner_triangle

# Area of one shaded region
area_one_region = area_segment + area_hat_triangle
# A simpler form is: Area(one region) = area_sector + area_polygon - 2 * area_inner_triangle
# Or Area(one region) = 1156*arctan(0.5) - 1734/5

# Step 3: Calculate the total area of both shaded regions.
total_area = 2 * area_one_region

print("The final calculation for the total shaded area is:")
print(f"Total Area = 2 * (Area of Sector - Area of Inner Triangle + Area of Hat Triangle)")
print(f"Total Area = 2 * (({radius**2} * arctan({side_length}/{radius})) - ({0.5*radius**2} * sin(2*arctan({side_length}/{radius}))) + (2*0.5*{radius}*{side_length} - {0.5*radius**2}*sin(2*arctan({side_length}/{radius}))))")
print(f"Total Area = 2 * (({radius**2} * math.atan({side_length / radius})) - {area_inner_triangle:.1f} + {area_hat_triangle:.1f})")
print(f"Total Area = {total_area:.4f}")
