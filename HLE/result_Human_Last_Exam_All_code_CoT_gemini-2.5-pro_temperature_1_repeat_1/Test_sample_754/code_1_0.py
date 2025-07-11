import math

# Step 1: Define constants and initial calculations
distance_ly = 10.0
speed_c = 0.02
travel_time_years = distance_ly / speed_c

initial_energy_per_ball_MJ = 100.0
leak_rate_per_year = 0.001
required_energy_MJ = 1000.0

box_dims = [12.0, 11.0, 11.0]
ball_radius = 2.0
ball_diameter = 2 * ball_radius
ball_diameter_sq = ball_diameter**2
grid_step = 0.5

# Calculate energy for a single leaking ball after the journey
energy_leaking_ball_MJ = initial_energy_per_ball_MJ * (1 - leak_rate_per_year) ** travel_time_years
energy_internal_ball_MJ = initial_energy_per_ball_MJ

print("--- Mission Analysis ---")
print(f"Travel Time: {travel_time_years:.0f} years")
print(f"Energy of a Non-Leaking Ball on Arrival: {energy_internal_ball_MJ:.2f} MJ")
print(f"Energy of a Leaking Ball on Arrival: {energy_leaking_ball_MJ:.2f} MJ")
print("-" * 30)

# Step 2: Calculate maximal number of balls using a greedy packing algorithm
print("Calculating maximum number of energy balls...")

# Define the valid space for the ball's center
center_min = [ball_radius, ball_radius, ball_radius]
center_max = [box_dims[0] - ball_radius, box_dims[1] - ball_radius, box_dims[2] - ball_radius]

# Helper function to generate a range of floats
def float_range(start, stop, step):
    while start < stop + (step / 2): # Add tolerance for float comparison
        yield round(start, 2) # Round to handle precision issues
        start += step

x_coords = list(float_range(center_min[0], center_max[0], grid_step))
y_coords = list(float_range(center_min[1], center_max[1], grid_step))
z_coords = list(float_range(center_min[2], center_max[2], grid_step))

placed_ball_centers = []

# Iterate through all possible center positions and place a ball if possible
for z in z_coords:
    for y in y_coords:
        for x in x_coords:
            potential_center = [x, y, z]
            is_valid_position = True
            for placed_center in placed_ball_centers:
                dist_sq = (potential_center[0] - placed_center[0])**2 + \
                          (potential_center[1] - placed_center[1])**2 + \
                          (potential_center[2] - placed_center[2])**2
                if dist_sq < ball_diameter_sq:
                    is_valid_position = False
                    break
            
            if is_valid_position:
                placed_ball_centers.append(potential_center)

total_balls = len(placed_ball_centers)
print(f"Maximal balls packed: {total_balls}")
print("-" * 30)

# Step 3: Count internal and touching balls
num_internal_balls = 0
num_touching_balls = 0

for center in placed_ball_centers:
    is_touching = (
        math.isclose(center[0], center_min[0]) or math.isclose(center[0], center_max[0]) or
        math.isclose(center[1], center_min[1]) or math.isclose(center[1], center_max[1]) or
        math.isclose(center[2], center_min[2]) or math.isclose(center[2], center_max[2])
    )
    if is_touching:
        num_touching_balls += 1
    else:
        num_internal_balls += 1

print(f"Number of internal (non-leaking) balls: {num_internal_balls}")
print(f"Number of touching (leaking) balls: {num_touching_balls}")
print("-" * 30)

# Step 4: Calculate total final energy
print("Calculating total available energy...")
total_energy = (num_internal_balls * energy_internal_ball_MJ) + (num_touching_balls * energy_leaking_ball_MJ)

# Print the final equation with all numbers
print("Final Energy Equation:")
print(f"({num_internal_balls} internal balls * {energy_internal_ball_MJ:.2f} MJ) + ({num_touching_balls} touching balls * {energy_leaking_ball_MJ:.2f} MJ) = {total_energy:.2f} MJ")
print("-" * 30)

# Step 5: Final conclusion
print("--- Final Conclusion ---")
is_sufficient = total_energy >= required_energy_MJ
final_answer = total_balls if is_sufficient else 0

print(f"Available Energy: {total_energy:.2f} MJ")
print(f"Required Energy: {required_energy_MJ:.2f} MJ")

if is_sufficient:
    print("Energy is SUFFICIENT for operations.")
else:
    print("Energy is NOT SUFFICIENT for operations.")

print(f"<<<{final_answer}>>>")