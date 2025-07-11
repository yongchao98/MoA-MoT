import math

# This script calculates the total energy for a proposed container design
# based on a heuristic optimization strategy.

# 1. Define the parameters of the problem
# Energy Balls
energy_large_ball = 20  # MJ for a 2-cm radius ball
energy_small_ball = 1   # MJ for a 1-cm radius ball
diameter_large_ball = 4.0 # cm
diameter_small_ball = 2.0 # cm

# Container Constraints
max_surface_area = 1050.0

# 2. Define the proposed container design
# Based on my analysis, a box is better suited for grid-packing.
# A 12x12x14 box is chosen to maximize large balls and then add small ones.
container_shape = "box"
L, W, H = 12.0, 12.0, 14.0

# 3. Calculate the number of packed balls
# Number of large balls is determined by how many fit in a simple grid
# within the container dimensions.
num_large_balls_L = math.floor(L / diameter_large_ball)
num_large_balls_W = math.floor(W / diameter_large_ball)
num_large_balls_H = math.floor(H / diameter_large_ball)
num_large_balls = num_large_balls_L * num_large_balls_W * num_large_balls_H

# The large balls occupy a 12x12x12 cm space. The remaining space is a
# 12x12x2 cm slab on top, which we fill with small balls.
remaining_height = H - (num_large_balls_H * diameter_large_ball)
num_small_balls_L = math.floor(L / diameter_small_ball)
num_small_balls_W = math.floor(W / diameter_small_ball)
num_small_balls_H = math.floor(remaining_height / diameter_small_ball)
num_small_balls = num_small_balls_L * num_small_balls_W * num_small_balls_H

# 4. Calculate total energy and verify surface area
total_energy = (num_large_balls * energy_large_ball) + (num_small_balls * energy_small_ball)
surface_area = 2 * (L*W + L*H + W*H)

# 5. Print the results of the analysis
print("--- Optimal Container Design Analysis ---")
print(f"Chosen Container: {container_shape.capitalize()} {int(L)}x{int(W)}x{int(H)} cm")
print(f"Calculated Surface Area: {surface_area:.2f} cm^2 (Constraint: <= {max_surface_area} cm^2)")
print(f"Surface area constraint is met: {surface_area <= max_surface_area}")
print("-" * 20)
print("Ball Packing Details:")
print(f"Number of 1-cm radius balls (a): {num_small_balls}")
print(f"Number of 2-cm radius balls (b): {num_large_balls}")
print("-" * 20)
print("Total Potential Energy Calculation:")
# Fulfilling the requirement to output each number in the final equation
print(f"Equation: ({num_large_balls} * {energy_large_ball} MJ) + ({num_small_balls} * {energy_small_ball} MJ)")
print(f"Total Energy: {total_energy} MJ")
print("-" * 20)

# The final answer in the required format [C]a;b
# C = container description, a = number of 1-cm balls, b = number of 2-cm balls
final_answer_format = f"[{container_shape} {int(L)}x{int(W)}x{int(H)}]{num_small_balls};{num_large_balls}"
print(f"Final Answer String: {final_answer_format}")
