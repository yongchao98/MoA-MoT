# Here is the Python code to calculate the number of turns.
# We analyze each move in the sequence to determine its contribution to the total count of right, left, and back presentations.

print("To solve this, we will analyze the turns in each dance move and sum them up.")
print("The counts are for turning the Right side / Left side / Back to the audience.")
print("-" * 20)

# Move Analysis:

# 1. Swivel Step Half Turn (16 counts): A 180-degree turn to the back, and a 180-degree turn to return to the front.
# This sequence presents the right side once, the left side once, and the back once.
swivel_half_turn_r, swivel_half_turn_l, swivel_half_turn_b = 1, 1, 1
print("Swivel Step Half Turn contributes: 1 right, 1 left, 1 back")

# 2. Sunanda (once): This move consists of two full 360-degree turns (e.g., one clockwise, one counter-clockwise).
# Each full turn presents the right side, back, and left side. Therefore, two turns contribute two of each.
sunanda_r, sunanda_l, sunanda_b = 2, 2, 2
print("Sunanda contributes: 2 right, 2 left, 2 back")

# 3. Barrel Turn (once): A fast 360-degree spin.
# The right and left sides are clearly presented. The back-facing moment is considered a transient part of the spin and not a full presentation.
barrel_turn_r, barrel_turn_l, barrel_turn_b = 1, 1, 0
print("Barrel Turn contributes: 1 right, 1 left, 0 back")
print("Other moves (Swivel Step, Balancing Step, Figure 8) involve no turns and contribute 0.")
print("-" * 20)
print("Calculating the total turns:")

# Sum the contributions from all turning moves.
# The moves with 0 turns are omitted from the sum for clarity.

# Calculate total right turns
total_right = swivel_half_turn_r + sunanda_r + barrel_turn_r
print(f"Total Right turns = {swivel_half_turn_r} + {sunanda_r} + {barrel_turn_r} = {total_right}")

# Calculate total left turns
total_left = swivel_half_turn_l + sunanda_l + barrel_turn_l
print(f"Total Left turns = {swivel_half_turn_l} + {sunanda_l} + {barrel_turn_l} = {total_left}")

# Calculate total back turns
total_back = swivel_half_turn_b + sunanda_b + barrel_turn_b
print(f"Total Back turns = {swivel_half_turn_b} + {sunanda_b} + {barrel_turn_b} = {total_back}")
print("-" * 20)

# Print the final result
print(f"The final count is {total_right} times turning the right side, {total_left} times turning the left side, and {total_back} times turning the back to the audience.")