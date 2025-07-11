import math

# Start of the Python script for the Interstellar problem.

# --- Problem Constants ---
DISTANCE_LY = 10.0  # light-years
SPEED_FRAC_C = 0.02 # fraction of speed of light
ENERGY_PER_BALL_MJ = 100.0
LEAK_RATE_PER_YEAR = 0.001 # 0.1%
REQUIRED_ENERGY_MJ = 1000.0
BALL_RADIUS_CM = 2.0
BOX_DIMS_CM = (12, 11, 11)

# --- Step 1: Calculate travel time ---
# The journey time is the distance divided by the speed.
travel_time_years = DISTANCE_LY / SPEED_FRAC_C
print("1. Calculating Travel Time")
print(f"The journey will take {DISTANCE_LY} light-years / {SPEED_FRAC_C}c = {travel_time_years:.0f} years.\n")

# --- Step 2: Calculate maximum number of balls ---
# The space available for ball centers is determined by the box dimensions minus the ball radius on each side.
# The number of balls that can be placed along an axis is based on the available space and the ball diameter.
ball_diameter = 2 * BALL_RADIUS_CM
Lx, Ly, Lz = BOX_DIMS_CM

# Calculate how many balls fit along each dimension in a grid packing
nx = math.floor((Lx - ball_diameter) / ball_diameter) + 1
ny = math.floor((Ly - ball_diameter) / ball_diameter) + 1
nz = math.floor((Lz - ball_diameter) / ball_diameter) + 1
max_balls = nx * ny * nz

print("2. Calculating Container Capacity")
print(f"Container dimensions: {Lx}x{Ly}x{Lz} cm. Ball diameter: {ball_diameter} cm.")
print(f"Max balls along each axis (x, y, z): {nx}, {ny}, {nz}")
print(f"The maximum number of balls that can be stored is {nx} * {ny} * {nz} = {max_balls}.\n")

# --- Step 3: Classify balls as touching or inner ---
# Inner balls are those not on the outer layer of the grid.
# A grid of nx*ny*nz balls has an inner core of (nx-2)*(ny-2)*(nz-2) balls.
inner_nx = max(0, nx - 2)
inner_ny = max(0, ny - 2)
inner_nz = max(0, nz - 2)
inner_balls = inner_nx * inner_ny * inner_nz
touching_balls = max_balls - inner_balls

print("3. Classifying the Energy Balls")
print(f"Number of inner balls (not leaking) = max(0, {nx}-2) * max(0, {ny}-2) * max(0, {nz}-2) = {inner_balls}")
print(f"Number of touching balls (leaking) = {max_balls} - {inner_balls} = {touching_balls}.\n")

# --- Step 4: Calculate total energy upon arrival ---
# Energy of a touching ball decays according to the compound decay formula.
energy_inner_final = ENERGY_PER_BALL_MJ
energy_touching_final = ENERGY_PER_BALL_MJ * (1 - LEAK_RATE_PER_YEAR) ** travel_time_years

# The total energy is the sum of energies from all inner and all touching balls.
total_energy_final = (inner_balls * energy_inner_final) + (touching_balls * energy_touching_final)

print("4. Calculating Final Energy")
print("Equation for a single touching ball's final energy:")
print(f"Final Energy = Initial Energy * (1 - Leak Rate) ^ Time")
print(f"Final Energy = {ENERGY_PER_BALL_MJ:.0f} MJ * (1 - {LEAK_RATE_PER_YEAR}) ^ {travel_time_years:.0f} = {energy_touching_final:.2f} MJ\n")

print("Equation for total final energy:")
print("Total Energy = (Inner Balls * Energy) + (Touching Balls * Final Energy)")
print(f"Total Energy = ({inner_balls} * {energy_inner_final:.2f}) + ({touching_balls} * {energy_touching_final:.2f}) = {total_energy_final:.2f} MJ\n")

# --- Step 5: Compare and decide the final answer ---
is_enough = total_energy_final >= REQUIRED_ENERGY_MJ

print("5. Final Conclusion")
print(f"The probe requires {REQUIRED_ENERGY_MJ} MJ upon arrival.")
print(f"The total energy available upon arrival is {total_energy_final:.2f} MJ.")

if is_enough:
    print("The energy is SUFFICIENT for the mission.")
    final_answer = max_balls
else:
    print("The energy is NOT SUFFICIENT for the mission.")
    final_answer = 0

print(f"\nSince the energy is insufficient, the answer is 0.")
print(f"<<<{final_answer}>>>")