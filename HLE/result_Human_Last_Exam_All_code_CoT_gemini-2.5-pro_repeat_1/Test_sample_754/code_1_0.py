import math

# Step 1: Define constants and calculate travel time
distance_ly = 10
speed_c = 0.02
travel_time_years = distance_ly / speed_c

initial_energy_per_ball_MJ = 100.0
leakage_rate_per_year = 0.001  # 0.1%
required_energy_MJ = 1000.0

box_dims = (12, 11, 11)
ball_radius = 2.0
ball_diameter = ball_radius * 2

# Step 2: Calculate the maximum number of balls that can be placed in a grid
nx = math.floor(box_dims[0] / ball_diameter)
ny = math.floor(box_dims[1] / ball_diameter)
nz = math.floor(box_dims[2] / ball_diameter)
max_balls = nx * ny * nz

# Step 3: Determine number of inner and surface balls
# A ball is an "inner" ball if it doesn't touch any of the 6 faces.
# With a simple grid packing starting at the edge, centers along the x-axis (length 12) are at 2, 6, 10.
# The layer at x=2 touches the x=0 face. The layer at x=10 touches the x=12 face. The layer at x=6 is an inner layer.
# Centers along the y-axis (length 11) are at 2, 6.
# The layer at y=2 touches the y=0 face. The layer at y=6 (ball extends to y=8) does NOT touch the y=11 face.
# Same for the z-axis (length 11).
num_inner_layers_x = nx - 2 if nx > 2 else 0
num_inner_layers_y = ny - 1 if ny > 1 else 0
num_inner_layers_z = nz - 1 if nz > 1 else 0

num_inner_balls = num_inner_layers_x * num_inner_layers_y * num_inner_layers_z
num_surface_balls = max_balls - num_inner_balls

# Step 4: Calculate the final energy of a single leaking ball
final_energy_leaking_ball = initial_energy_per_ball_MJ * math.pow(1 - leakage_rate_per_year, travel_time_years)

# Step 5: Calculate the total final energy
total_energy = (num_inner_balls * initial_energy_per_ball_MJ) + (num_surface_balls * final_energy_leaking_ball)

# Step 6: Compare, print results, and conclude
print(f"The probe will travel for {int(travel_time_years)} years.")
print(f"The container can hold a maximum of {max_balls} balls, arranged in a {nx}x{ny}x{nz} grid.")
print(f"Based on this packing, there are {num_inner_balls} inner ball(s) and {num_surface_balls} surface ball(s).")
print("\nThe final energy is calculated as follows:")
print(f"Total Energy = (Number of inner balls * Initial Energy) + (Number of surface balls * Final Energy of a leaking ball)")
print(f"Total Energy = ({num_inner_balls} * {initial_energy_per_ball_MJ}) + ({num_surface_balls} * {initial_energy_per_ball_MJ} * (1 - {leakage_rate_per_year})^{int(travel_time_years)})")
print(f"Total Energy = {num_inner_balls * initial_energy_per_ball_MJ} + {num_surface_balls} * {final_energy_leaking_ball:.3f}")
print(f"Total Energy = {total_energy:.3f} MJ")

print(f"\nThe required energy for operations is {required_energy_MJ} MJ.")

if total_energy >= required_energy_MJ:
    final_answer = max_balls
    print(f"The total available energy ({total_energy:.3f} MJ) is sufficient.")
else:
    final_answer = 0
    print(f"The total available energy ({total_energy:.3f} MJ) is NOT sufficient.")

print(f"\nFinal Answer: {final_answer}")

print("<<<0>>>")