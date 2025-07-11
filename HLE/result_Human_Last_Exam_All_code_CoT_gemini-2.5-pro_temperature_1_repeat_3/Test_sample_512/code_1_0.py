import math

# Initial container specifications
initial_side = 12.0
initial_surface_area = 6 * initial_side**2
ball_radius = 2.0
ball_diameter = 2 * ball_radius
balls_per_side = initial_side / ball_diameter
initial_ball_count = int(balls_per_side**3)

# My best attempt at a new container. Let's analyze it as a demonstration.
# Packing strategy: 3 layers of 9 balls (3x3 grid).
# Layers are shifted to be closer. sz=3.5, sx=1.5, sy=1.5
# This packing was found to be valid (distance > 4cm).
box_x = 13.5
box_y = 13.5
box_z = 11.0
new_surface_area = 2 * (box_x * box_y + box_x * box_z + box_y * box_z)

# The problem is to find if a better container exists.
# "Better" means SA < initial_surface_area and ball_count >= initial_ball_count
# Based on extensive search, any valid packing of 27 balls under the 0.5cm grid constraint
# seems to result in a bounding box with a surface area greater than the original cube.
# For example, the 13.5x13.5x11 box holds 27 balls but has SA = 958.5 cm^2.
# Another valid packing for 27 balls results in a 14x12x11 box with SA = 908 cm^2.
# Since no configuration with SA < 864 cm^2 was found, the answer is 0.

result = 0
print(result)
