import math

# --- Problem Parameters ---
# Number of 2-cm radius balls (large)
num_large_balls = 27
# Number of 1-cm radius balls (small)
num_small_balls = 44

# Energy per ball in MJ
energy_large_ball = 20
energy_small_ball = 1

# Container description
container_description = "box 12x12x15"

# --- Calculation ---
# Total energy from large balls
total_energy_large = num_large_balls * energy_large_ball
# Total energy from small balls
total_energy_small = num_small_balls * energy_small_ball
# Total energy
total_energy = total_energy_large + total_energy_small

# --- Output ---
print(f"The best container found is a {container_description}.")
print(f"It can hold {num_large_balls} large balls and {num_small_balls} small balls.")
print("The total energy calculation is:")
print(f"{num_large_balls} * {energy_large_ball} MJ + {num_small_balls} * {energy_small_ball} MJ = {total_energy} MJ")
