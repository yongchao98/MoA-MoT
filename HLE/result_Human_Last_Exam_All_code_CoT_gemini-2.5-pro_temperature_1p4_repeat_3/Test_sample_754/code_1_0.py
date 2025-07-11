import math

# Step 1: Define constants and calculate travel time
distance_ly = 10.0  # light-years
speed_c = 0.02  # speed as a fraction of c
travel_time_years = distance_ly / speed_c

# Step 2: Calculate the final energy of a single leaking ball
initial_energy_mj = 100.0
leak_rate_per_year = 0.001  # 0.1%
leak_multiplier = (1 - leak_rate_per_year) ** travel_time_years
energy_surface_ball_final = initial_energy_mj * leak_multiplier
energy_internal_ball_final = initial_energy_mj

# Step 3: Define the optimal 16-ball packing configuration found
# This packing consists of 2 layers, each with 8 balls.
# Layer z-coordinates
z_coords = [2.0, 6.0]
# Layer (x,y) coordinates for the 8 balls
layer_xy_coords = [
    (2.0, 2.0), (6.0, 2.0), (10.0, 2.0),  # Row 1
    (4.0, 5.5), (8.0, 5.5),              # Row 2
    (2.0, 9.0), (6.0, 9.0), (10.0, 9.0)   # Row 3
]

all_ball_centers = []
for z in z_coords:
    for xy in layer_xy_coords:
        all_ball_centers.append((xy[0], xy[1], z))

# Step 4: Count internal and surface balls for this configuration
box_dims = (12.0, 11.0, 11.0)
ball_radius = 2.0
n_internal = 0
n_surface = 0

for center in all_ball_centers:
    x, y, z = center
    # A ball is internal if it doesn't touch any wall
    # Condition for touching: center_coord <= radius or center_coord >= dimension - radius
    is_surface = (
        x <= ball_radius or x >= box_dims[0] - ball_radius or
        y <= ball_radius or y >= box_dims[1] - ball_radius or
        z <= ball_radius or z >= box_dims[2] - ball_radius
    )
    if is_surface:
        n_surface += 1
    else:
        n_internal += 1

# Step 5: Calculate total final energy
total_energy = (n_internal * energy_internal_ball_final) + \
               (n_surface * energy_surface_ball_final)

# Step 6: Check if energy is sufficient and output the result
energy_requirement_mj = 1000.0

print(f"Journey time: {travel_time_years} years")
print(f"Final energy of one internal ball: {energy_internal_ball_final:.2f} MJ")
print(f"Final energy of one surface ball: {energy_surface_ball_final:.2f} MJ")
print("-" * 30)
print(f"Proposed packing configuration: {len(all_ball_centers)} balls")
print(f"Number of internal balls (N_internal): {n_internal}")
print(f"Number of surface balls (N_surface): {n_surface}")
print("-" * 30)
print("Final energy equation:")
print(f"{n_internal} * {energy_internal_ball_final:.2f} MJ + {n_surface} * {energy_surface_ball_final:.2f} MJ = {total_energy:.2f} MJ")

if total_energy >= energy_requirement_mj:
    print(f"\nThe total energy ({total_energy:.2f} MJ) is sufficient (>= {energy_requirement_mj} MJ).")
    answer = n_internal + n_surface
else:
    print(f"\nThe total energy ({total_energy:.2f} MJ) is insufficient (< {energy_requirement_mj} MJ).")
    answer = 0
    
print("\nIs the container enough? If yes, what is the maximal number of balls? (0 if no)")
print(f"<<<{answer}>>>")