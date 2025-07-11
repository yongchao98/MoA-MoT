import math

# --- Problem Parameters ---
distance_ly = 10.0
speed_c = 0.02
initial_energy_per_ball = 100.0  # MJ
leak_rate_per_year = 0.001  # 0.1%
required_energy = 1000.0  # MJ
ball_radius = 2.0  # cm
ball_diameter = 4.0 # cm
box_dims = (12.0, 11.0, 11.0)

# --- Step 1: Calculate travel time ---
travel_time_years = distance_ly / speed_c
print(f"Step 1: Calculate travel time")
print(f"The probe travels 10 light-years at 0.02c.")
print(f"Travel time = {distance_ly:.0f} light-years / {speed_c}c = {travel_time_years:.0f} years")
print("-" * 40)

# --- Step 2: Calculate maximum number of balls that can fit ---
# We place the balls in a simple cubic grid.
# For a ball with radius 2cm to fit, its center (x,y,z) must be within:
# x: [2, 10], y: [2, 9], z: [2, 9]
# With a 4cm diameter, we can place centers at:
# x-coords: 2, 6, 10 (3 balls)
# y-coords: 2, 6 (2 balls)
# z-coords: 2, 6 (2 balls)
num_x = 3
num_y = 2
num_z = 2
total_balls = num_x * num_y * num_z
print(f"Step 2: Calculate the maximum number of balls")
print(f"Container: {box_dims[0]}x{box_dims[1]}x{box_dims[2]} cm. Ball diameter: {ball_diameter} cm.")
print(f"Max balls that can be packed = {num_x} * {num_y} * {num_z} = {total_balls}")
print("-" * 40)

# --- Step 3: Identify leaking and non-leaking balls ---
# A ball leaks if it touches the container surface.
# A ball with radius 2cm is touching if its center's x-coord is 2 or 10, y-coord is 2, or z-coord is 2.
# The only center coordinate combination that doesn't meet these conditions is (6, 6, 6).
num_non_leaking = 1
num_leaking = total_balls - num_non_leaking
print(f"Step 3: Identify leaking vs. non-leaking balls")
print(f"Out of {total_balls} balls, only the one with center (6.0, 6.0, 6.0) does not touch any surface.")
print(f"Number of leaking balls: {num_leaking}")
print(f"Number of non-leaking balls: {num_non_leaking}")
print("-" * 40)

# --- Step 4: Calculate the final energy ---
# Energy of a leaking ball after the journey
energy_leaking_ball = initial_energy_per_ball * (1 - leak_rate_per_year) ** travel_time_years
# Energy of a non-leaking ball is unchanged
energy_non_leaking_ball = initial_energy_per_ball
total_energy = (num_leaking * energy_leaking_ball) + (num_non_leaking * energy_non_leaking_ball)

print(f"Step 4: Calculate total energy upon arrival")
print(f"Energy per leaking ball = {initial_energy_per_ball:.0f} * (1 - {leak_rate_per_year}) ^ {travel_time_years:.0f} = {energy_leaking_ball:.2f} MJ")
print(f"Energy per non-leaking ball = {energy_non_leaking_ball:.0f} MJ")
print("\nFinal Energy Calculation:")
print(f"Total Energy = ({num_leaking} * {energy_leaking_ball:.2f}) + ({num_non_leaking} * {energy_non_leaking_ball:.0f})")
print(f"Total Energy = {num_leaking * energy_leaking_ball:.2f} MJ + {num_non_leaking * energy_non_leaking_ball:.0f} MJ = {total_energy:.2f} MJ")
print("-" * 40)

# --- Step 5: Compare with requirement and give the answer ---
print(f"Step 5: Final conclusion")
print(f"Required energy for operations: {required_energy:.0f} MJ.")
print(f"Available energy upon arrival: {total_energy:.2f} MJ.")

if total_energy >= required_energy:
    print("The available energy is sufficient.")
    answer = total_balls
else:
    print("The available energy is NOT sufficient.")
    answer = 0

print(f"<<<{answer}>>>")