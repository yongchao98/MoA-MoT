import math

# Define problem constants
ENERGY_REQ = 1000  # MJ
ENERGY_PER_BALL = 25  # MJ
COST_PER_BALL = 1000  # USD
MATERIAL_COST_PER_CM2 = 200  # USD/cm^2
BALL_RADIUS = 2  # cm
BALL_DIAMETER = 4 # cm
PRECISION = 0.5  # cm

# Optimal design parameters determined from analysis
# This configuration uses 42 balls, which exceeds the minimum requirement of 40,
# but results in a lower total cost due to more efficient packing.
balls_per_layer = 7
num_layers = 6

# --- Calculations for the Optimal Cylindrical Design ---

# 1. Calculate the total number of balls and the cost of the balls
total_balls = balls_per_layer * num_layers
cost_of_balls = total_balls * COST_PER_BALL

# 2. Calculate the dimensions of the cylinder container
# For a layer of 7 balls (1 in center, 6 around), the required radius is 3 * ball_radius.
actual_radius = 3 * BALL_RADIUS
# The container radius must be a multiple of the precision.
container_radius = math.ceil(actual_radius / PRECISION) * PRECISION

# For height, assume densest stacking (like HCP), where the height added by each layer
# after the first is diameter * sqrt(2/3).
actual_height = (num_layers - 1) * BALL_DIAMETER * math.sqrt(2.0/3.0) + BALL_DIAMETER
# The container height must be a multiple of the precision.
container_height = math.ceil(actual_height / PRECISION) * PRECISION

# 3. Calculate the surface area of the cylinder
surface_area_top_bottom = 2 * math.pi * container_radius**2
surface_area_side = 2 * math.pi * container_radius * container_height
total_surface_area = surface_area_top_bottom + surface_area_side

# 4. Calculate the cost of the container material
cost_of_material = total_surface_area * MATERIAL_COST_PER_CM2

# 5. Calculate the total cost
total_cost = cost_of_balls + cost_of_material

# --- Output the results ---
print("Optimal Design: Cylinder")
print(f"Packing Arrangement: {num_layers} layers of {balls_per_layer} balls each")
print(f"Total Energy Balls: {total_balls}")
print(f"Container Radius: {container_radius} cm")
print(f"Container Height: {container_height} cm")
print(f"Container Surface Area: {total_surface_area:.2f} cm^2")
print("\n--- Cost Breakdown ---")
print(f"Cost of Energy Balls: ${int(cost_of_balls)}")
print(f"Cost of Container Material: ${cost_of_material:.2f}")

print("\n--- Final Total Cost Calculation ---")
# Output each number in the final equation as requested
print(f"C = ${cost_of_balls:.2f} (Balls) + ${cost_of_material:.2f} (Material)")
print(f"C = ${total_cost:.2f}")

# Final answer format
print(f"\n<<<{total_cost:.2f}>>>")