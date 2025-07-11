import math

# --- Problem Parameters ---
# Energy Ball properties
ball_radius_cm = 2.0
ball_diameter_cm = ball_radius_cm * 2
energy_per_ball_MJ = 30
cost_per_ball_usd = 1000

# Container properties
max_surface_area_cm2 = 1000.0
container_material_cost_per_cm2 = 200

# Mission requirements
min_total_energy_MJ = 1000.0

# --- Step 1: Calculate Minimum Number of Balls ---
# We need at least 1000 MJ of energy.
num_balls = math.ceil(min_total_energy_MJ / energy_per_ball_MJ)

# --- Step 2: Determine Optimal Container Design ---
# Based on the analysis, a box is not feasible as its surface area exceeds 1000 cm^2.
# The best feasible design is a cylinder capable of holding 5 layers of 7 balls each.

# Optimal cylinder dimensions
# A layer of 7 balls (1 center, 6 around) of radius 2cm fits in a circle of radius 6cm.
cyl_radius_cm = 6.0
# To hold 34 balls (ceil(34/7)=5), we need 5 layers. Each layer is 4cm high.
cyl_height_cm = 5 * ball_diameter_cm

# --- Step 3: Calculate Costs for the Chosen Design ---
# Calculate the container's surface area
surface_area_top_bottom = 2 * math.pi * cyl_radius_cm**2
surface_area_side = 2 * math.pi * cyl_radius_cm * cyl_height_cm
total_surface_area = surface_area_top_bottom + surface_area_side

# Calculate the cost of the container
container_cost = total_surface_area * container_material_cost_per_cm2

# Calculate the cost of the energy balls
balls_cost = num_balls * cost_per_ball_usd

# Calculate the total cost
total_cost = container_cost + balls_cost

# --- Step 4: Output the Result ---
# Print the breakdown of the final calculation
print("--- Optimal Design: Cylinder Container ---")
print(f"Number of Energy Balls: {int(num_balls)}")
print(f"Container Dimensions: Radius = {cyl_radius_cm} cm, Height = {cyl_height_cm} cm")
print(f"Container Surface Area: {total_surface_area:.2f} cm^2\n")

print("--- Cost Calculation ---")
print("Total Cost = (Container Surface Area * Cost per cm^2) + (Number of Balls * Cost per Ball)")
print(f"Total Cost = ({total_surface_area:.4f} * {container_material_cost_per_cm2}) + ({int(num_balls)} * {cost_per_ball_usd})")
print(f"Total Cost = {container_cost:.2f} + {balls_cost:.2f}")
print(f"Total Cost = {total_cost:.2f} USD")

# Final answer in the required format
print(f"\n<<<{total_cost}>>>")