import math

# Step 1: Define constants and calculate mission duration
distance_ly = 10  # light-years
speed_c = 0.02    # proportion of speed of light
travel_time_years = distance_ly / speed_c

# Step 2: Calculate final energy per ball
initial_energy_MJ = 100
leakage_rate_per_year = 0.001 # 0.1%
total_leakage = leakage_rate_per_year * travel_time_years

# Energy of a ball that touches the container walls
energy_touching_ball = initial_energy_MJ * (1 - total_leakage)
# Energy of an inner ball (no leakage)
energy_inner_ball = initial_energy_MJ

# Step 3: Determine maximum ball capacity
box_dims = [12, 11, 11] # cm
ball_radius = 2 # cm
ball_diameter = ball_radius * 2

# Calculate how many balls fit along each dimension using simple grid packing
# The number of balls is floor((dimension - diameter) / diameter) + 1, which simplifies to floor(dimension/diameter).
# However, the space available for centers is (dim - 2*r), so number of balls is floor((dim - 2*r)/d) + 1
balls_x = math.floor((box_dims[0] - 2 * ball_radius) / ball_diameter) + 1
balls_y = math.floor((box_dims[1] - 2 * ball_radius) / ball_diameter) + 1
balls_z = math.floor((box_dims[2] - 2 * ball_radius) / ball_diameter) + 1
total_balls = balls_x * balls_y * balls_z

# Step 4: Count touching vs. non-touching balls
# A ball is "inner" if it is not in the first or last layer along any axis.
# The packing grid starts with centers at x=2, y=2, z=2.
# X-centers: 2, 6, 10. Inner ball is at x=6.
# Y-centers: 2, 6. The ball at y=2 touches the wall at y=0. The ball at y=6 extends from y=4 to y=8, not touching the wall at y=11.
# Z-centers: 2, 6. The ball at z=2 touches. The ball at z=6 does not.
# An "inner" ball must not touch any wall. This is only true for the ball at center (6, 6, 6).
num_inner_balls_x = 1 # The one at center x=6
num_inner_balls_y = 1 # The one at center y=6
num_inner_balls_z = 1 # The one at center z=6
num_inner_balls = num_inner_balls_x * num_inner_balls_y * num_inner_balls_z

num_touching_balls = total_balls - num_inner_balls

# Step 5: Calculate total available energy
total_final_energy = (num_touching_balls * energy_touching_ball) + (num_inner_balls * energy_inner_ball)

# Step 6: Compare and provide the answer
required_energy_MJ = 1000

print(f"Pioneer's mission to Pandora will take {int(travel_time_years)} years.")
print(f"The container can hold a maximum of {total_balls} energy balls ({balls_x}x{balls_y}x{balls_z}).")
print(f"After {int(travel_time_years)} years:")
print(f"- Each of the {num_touching_balls} touching balls will have {energy_touching_ball:.1f} MJ.")
print(f"- Each of the {num_inner_balls} inner ball will have {energy_inner_ball:.1f} MJ.")
print("The final energy calculation is:")
print(f"({num_touching_balls} * {energy_touching_ball:.1f}) + ({num_inner_balls} * {energy_inner_ball:.1f}) = {total_final_energy:.1f} MJ")

print("\n--- Final Answer ---")
# If the total energy is sufficient, the answer is the number of balls. Otherwise, it's 0.
if total_final_energy >= required_energy_MJ:
    result = total_balls
else:
    result = 0

print(f"The total available energy ({total_final_energy:.1f} MJ) is less than the required {required_energy_MJ} MJ.")
print(f"Therefore, the probe will not have enough energy for operations.")
print(f"The number of balls to be sent is {result}.")

# Final answer format
print(f"\n<<<{result}>>>")