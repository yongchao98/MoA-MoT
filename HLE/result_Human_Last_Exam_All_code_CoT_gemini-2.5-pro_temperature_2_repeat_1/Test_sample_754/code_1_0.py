import math

# Step 1: Define constants and calculate travel time
distance_ly = 10.0
speed_fraction_c = 0.02
ball_radius_cm = 2.0
ball_diameter_cm = ball_radius_cm * 2
initial_energy_mj = 100.0
leak_rate_per_year = 0.001  # 0.1%
container_dims_cm = (12.0, 11.0, 11.0)
required_energy_mj = 1000.0

# Time = Distance / Speed
travel_time_years = distance_ly / speed_fraction_c

# Step 2: Determine the maximum number of balls that can be placed
# The effective length for placing balls is Dimension - 2 * Radius
# We use floor division to find how many diameters fit, and add 1 for the first ball.
num_balls_x = math.floor((container_dims_cm[0] - 2 * ball_radius_cm) / ball_diameter_cm) + 1
num_balls_y = math.floor((container_dims_cm[1] - 2 * ball_radius_cm) / ball_diameter_cm) + 1
num_balls_z = math.floor((container_dims_cm[2] - 2 * ball_radius_cm) / ball_diameter_cm) + 1
max_total_balls = num_balls_x * num_balls_y * num_balls_z

# Step 3: Identify which balls are leaking
# A ball is leaking if it touches any of the 6 container walls.
# We find the number of non-leaking (interior) balls and subtract from the total.
# An interior ball's center (x,y,z) must not be at the minimum (radius) or maximum (dim-radius) possible position.

# Center positions for balls along each axis
center_positions_x = [ball_radius_cm + i * ball_diameter_cm for i in range(num_balls_x)]
center_positions_y = [ball_radius_cm + i * ball_diameter_cm for i in range(num_balls_y)]
center_positions_z = [ball_radius_cm + i * ball_diameter_cm for i in range(num_balls_z)]

# Count the number of non-leaking balls
num_non_leaking_balls = 0
for x in center_positions_x:
    for y in center_positions_y:
        for z in center_positions_z:
            # Check if the ball at center (x,y,z) is touching a wall
            is_leaking = (x == ball_radius_cm or x == container_dims_cm[0] - ball_radius_cm or
                          y == ball_radius_cm or y == container_dims_cm[1] - ball_radius_cm or
                          z == ball_radius_cm or z == container_dims_cm[2] - ball_radius_cm)
            if not is_leaking:
                num_non_leaking_balls += 1

num_leaking_balls = max_total_balls - num_non_leaking_balls

# Step 4 & 5: Calculate the final total energy
energy_per_leaking_ball = initial_energy_mj * ((1 - leak_rate_per_year) ** travel_time_years)
total_final_energy = (num_leaking_balls * energy_per_leaking_ball) + (num_non_leaking_balls * initial_energy_mj)

# Print the breakdown of the calculation
print("Problem Analysis:")
print(f" - Travel time: {int(travel_time_years)} years")
print(f" - Container can hold a maximum of {max_total_balls} balls ({num_balls_x}x{num_balls_y}x{num_balls_z}).")
print(f" - Balls touching container surface (leaking): {num_leaking_balls}")
print(f" - Balls not touching surface (non-leaking): {num_non_leaking_balls}")
print(f" - Energy per leaking ball on arrival: {energy_per_leaking_ball:.2f} MJ")
print("-" * 20)

# Print the final energy equation as requested
print("Final Energy Calculation:")
# Note: the equation shows the numbers used for the final calculation.
print(f"{num_leaking_balls} * {initial_energy_mj} * (1 - {leak_rate_per_year})^{int(travel_time_years)} + {num_non_leaking_balls} * {initial_energy_mj} = {total_final_energy:.2f} MJ")
print("-" * 20)

# Step 6: Make the final decision
print("Conclusion:")
print(f"The probe requires {required_energy_mj} MJ but will have {total_final_energy:.2f} MJ.")

if total_final_energy >= required_energy_mj:
    final_answer = max_total_balls
    print("The energy supply is sufficient.")
else:
    final_answer = 0
    print("The energy supply is NOT sufficient.")

# Output the final answer in the specified format
print(f"<<<{final_answer}>>>")