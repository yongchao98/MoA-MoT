import math

# Step 1: Define constants and calculate travel time
distance_ly = 10  # light-years
speed_c = 0.02    # fraction of speed of light
travel_time_years = distance_ly / speed_c

initial_energy_per_ball_MJ = 100
leak_rate_per_year = 0.001  # 0.1%
required_energy_MJ = 1000

# Box and ball dimensions
box_dims = (12, 11, 11)  # cm
ball_radius = 2          # cm
ball_diameter = 4        # cm

print(f"The probe will travel for {travel_time_years:.0f} years.")

# Step 2: Calculate energy remaining for a single leaking ball
energy_per_leaking_ball = initial_energy_per_ball_MJ * math.pow((1 - leak_rate_per_year), travel_time_years)
energy_per_non_leaking_ball = initial_energy_per_ball_MJ

print(f"Energy of a single non-leaking ball upon arrival: {energy_per_non_leaking_ball:.1f} MJ")
print(f"Energy of a single leaking ball upon arrival: {energy_per_leaking_ball:.3f} MJ")

# Step 3: Determine ball packing and count leaking/non-leaking balls
# We use a staggered packing arrangement to fit more balls than a simple cubic lattice.
# Layer 1 (z=2.0): A 3x2 grid of 6 balls
# Layer 2 (z=5.0): 2 balls placed in the hollows of layer 1
# Layer 3 (z=8.0): Another 3x2 grid of 6 balls
# This gives a total of 6 + 2 + 6 = 14 balls.

centers = [
    # Layer at z=2.0 (6 balls)
    (2.0, 2.0, 2.0), (6.0, 2.0, 2.0), (10.0, 2.0, 2.0),
    (2.0, 6.0, 2.0), (6.0, 6.0, 2.0), (10.0, 6.0, 2.0),
    # Layer at z=5.0 (2 balls)
    (4.0, 4.0, 5.0), (8.0, 4.0, 5.0),
    # Layer at z=8.0 (6 balls)
    (2.0, 2.0, 8.0), (6.0, 2.0, 8.0), (10.0, 2.0, 8.0),
    (2.0, 6.0, 8.0), (6.0, 6.0, 8.0), (10.0, 6.0, 8.0)
]

max_balls = len(centers)
print(f"\nA staggered packing configuration allows for a maximum of {max_balls} balls.")

# A ball leaks if its center is on the boundary of the allowed space for centers.
# Allowed center coordinates: cx in [2, 10], cy in [2, 9], cz in [2, 9].
# Boundary values are cx=2 or 10; cy=2 or 9; cz=2 or 9.

num_leaking = 0
num_non_leaking = 0
for cx, cy, cz in centers:
    is_leaking = False
    if cx == 2.0 or cx == 10.0:
        is_leaking = True
    if cy == 2.0 or cy == 9.0:
        is_leaking = True
    if cz == 2.0 or cz == 9.0:
        is_leaking = True
    
    if is_leaking:
        num_leaking += 1
    else:
        num_non_leaking += 1

print(f"Number of non-leaking balls: {num_non_leaking}")
print(f"Number of leaking balls: {num_leaking}")

# Step 4: Calculate total energy and check if it's sufficient
total_energy = (num_non_leaking * energy_per_non_leaking_ball) + \
               (num_leaking * energy_per_leaking_ball)

print("\nFinal energy calculation:")
print(f"Total Energy = (Non-leaking balls * Energy) + (Leaking balls * Energy)")
print(f"Total Energy = ({num_non_leaking} * {energy_per_non_leaking_ball:.1f}) + ({num_leaking} * {energy_per_leaking_ball:.3f})")
print(f"Total Energy = {num_non_leaking * energy_per_non_leaking_ball:.1f} + {num_leaking * energy_per_leaking_ball:.3f}")
print(f"Total Energy = {total_energy:.3f} MJ")

print(f"\nRequired energy for operations is {required_energy_MJ} MJ.")
if total_energy >= required_energy_MJ:
    print("The available energy is sufficient.")
    final_answer = max_balls
else:
    print("The available energy is NOT sufficient.")
    final_answer = 0

# Final Answer
print(f"The final answer is: {final_answer}")