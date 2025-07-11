import math

# Step 1: Define constants and calculate mission duration
distance_ly = 10  # light-years
speed_fraction_c = 0.02  # 0.02c
travel_time_years = distance_ly / speed_fraction_c

print(f"1. Mission Duration Calculation:")
print(f"   - Distance: {distance_ly} light-years")
print(f"   - Speed: {speed_fraction_c}c")
print(f"   - Travel Time = {distance_ly} / {speed_fraction_c} = {travel_time_years:.0f} years")
print("-" * 30)

# Step 2: Calculate energy per ball on arrival
initial_energy_MJ = 100.0
leak_rate_per_year = 0.001  # 0.1%
ball_radius_cm = 2.0
ball_diameter_cm = ball_radius_cm * 2

energy_loss_per_leaking_ball = initial_energy_MJ * leak_rate_per_year * travel_time_years
final_energy_leaking_ball = initial_energy_MJ - energy_loss_per_leaking_ball
final_energy_inner_ball = initial_energy_MJ

print(f"2. Energy Per Ball Calculation:")
print(f"   - Initial energy per ball: {initial_energy_MJ} MJ")
print(f"   - Energy loss for a leaking ball over {travel_time_years:.0f} years = {initial_energy_MJ} MJ * {leak_rate_per_year} * {travel_time_years:.0f} years = {energy_loss_per_leaking_ball:.1f} MJ")
print(f"   - Final energy of a leaking ball: {initial_energy_MJ} - {energy_loss_per_leaking_ball:.1f} = {final_energy_leaking_ball:.1f} MJ")
print(f"   - Final energy of an inner ball (no leaks): {final_energy_inner_ball:.1f} MJ")
print("-" * 30)

# Step 3: Calculate maximum ball capacity
container_dims_cm = (12, 11, 11)
L, W, H = container_dims_cm

# Number of balls along each dimension
# Formula: floor((Dimension - Diameter) / Diameter) + 1, equivalent to floor((Dim - 2*Rad) / (2*Rad)) + 1
# This calculates how many spheres of a given diameter can be placed starting from one edge.
nx = math.floor((L - 2 * ball_radius_cm) / ball_diameter_cm) + 1
ny = math.floor((W - 2 * ball_radius_cm) / ball_diameter_cm) + 1
nz = math.floor((H - 2 * ball_radius_cm) / ball_diameter_cm) + 1
total_balls = nx * ny * nz

print(f"3. Container Capacity Calculation:")
print(f"   - Container dimensions: {L}x{W}x{H} cm")
print(f"   - Ball diameter: {ball_diameter_cm} cm")
print(f"   - Balls along length ({L} cm): {nx}")
print(f"   - Balls along width ({W} cm): {ny}")
print(f"   - Balls along height ({H} cm): {nz}")
print(f"   - Total number of balls = {nx} * {ny} * {nz} = {total_balls}")
print("-" * 30)

# Step 4: Distinguish leaking vs. inner balls
# A ball is 'inner' if it's not on the outer layer of the 3D grid of balls.
# The number of inner balls is (nx-2) * (ny-2) * (nz-2). Must handle cases where ni-2 is negative.
inner_nx = max(0, nx - 2)
inner_ny = max(0, ny - 2)
inner_nz = max(0, nz - 2)
num_inner_balls = inner_nx * inner_ny * inner_nz
num_leaking_balls = total_balls - num_inner_balls

print(f"4. Leaking vs. Inner Balls Calculation:")
print(f"   - The balls form a {nx}x{ny}x{nz} grid.")
print(f"   - The inner core of this grid is {inner_nx}x{inner_ny}x{inner_nz}.")
print(f"   - Number of inner (non-leaking) balls = {num_inner_balls}")
print(f"   - Number of leaking balls = {total_balls} - {num_inner_balls} = {num_leaking_balls}")
print("-" * 30)

# Step 5: Calculate total final energy
total_final_energy = (num_inner_balls * final_energy_inner_ball) + (num_leaking_balls * final_energy_leaking_ball)

print("5. Total Final Energy Calculation:")
equation = f"   Total Energy = ({num_inner_balls} inner * {final_energy_inner_ball:.1f} MJ) + ({num_leaking_balls} leaking * {final_energy_leaking_ball:.1f} MJ)"
print(equation)
print(f"   Total Energy = {num_inner_balls * final_energy_inner_ball:.1f} MJ + {num_leaking_balls * final_energy_leaking_ball:.1f} MJ = {total_final_energy:.1f} MJ")
print("-" * 30)

# Step 6: Compare and conclude
required_energy_MJ = 1000.0
is_sufficient = total_final_energy >= required_energy_MJ
result = total_balls if is_sufficient else 0

print("6. Final Conclusion:")
print(f"   - Energy required upon arrival: {required_energy_MJ:.1f} MJ")
print(f"   - Energy available upon arrival: {total_final_energy:.1f} MJ")

if is_sufficient:
    print(f"   - The available energy is sufficient. The maximal number of balls is {result}.")
else:
    print(f"   - The available energy is NOT sufficient. The result is {result}.")

print("<<<0>>>")