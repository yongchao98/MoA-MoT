# Plan:
# 1. Define the number of right, left, and back turns for each dance move.
# 2. A key assumption is that "Swivel Step Half Turn (16 counts)" means the move is performed twice,
#    as a standard phrase for this move is 8 counts.
# 3. Sum the turns from all moves to get the total for each category.
# 4. Print the breakdown of the calculation and the final result.

# --- Step 1: Define turns for each move ---

# Swivel Step (8 counts): No body turns.
swivel_step_R, swivel_step_L, swivel_step_B = 0, 0, 0

# Swivel Step Half Turn (16 counts): Assumed to be two 8-count half-turns.
# Each half-turn presents the back to the audience.
swivel_step_half_turn_R, swivel_step_half_turn_L, swivel_step_half_turn_B = 0, 0, 2

# Sunanda (once): A turn away and a turn back. Presents right side, left side, and back once each.
sunanda_R, sunanda_L, sunanda_B = 1, 1, 1

# Balancing Step (once): No body turns.
balancing_step_R, balancing_step_L, balancing_step_B = 0, 0, 0

# Figure 8 (8 times): Standard hip move with no body turns.
figure_8_R, figure_8_L, figure_8_B = 0, 0, 0

# Barrel Turn (once): A full 360-degree spin presents each orientation once.
barrel_turn_R, barrel_turn_L, barrel_turn_B = 1, 1, 1

# --- Step 2: Calculate the totals ---

total_R = swivel_step_R + swivel_step_half_turn_R + sunanda_R + balancing_step_R + figure_8_R + barrel_turn_R
total_L = swivel_step_L + swivel_step_half_turn_L + sunanda_L + balancing_step_L + figure_8_L + barrel_turn_L
total_B = swivel_step_B + swivel_step_half_turn_B + sunanda_B + balancing_step_B + figure_8_B + barrel_turn_B

# --- Step 3: Print the results ---

print("Calculating the total number of turns based on the dance sequence analysis:")

# Print the equation for Right side turns
print(f"\nTotal Right Side Turns = (Swivel Step) {swivel_step_R} + (SSHT) {swivel_step_half_turn_R} + (Sunanda) {sunanda_R} + (Balancing) {balancing_step_R} + (Figure 8) {figure_8_R} + (Barrel) {barrel_turn_R} = {total_R}")

# Print the equation for Left side turns
print(f"Total Left Side Turns  = (Swivel Step) {swivel_step_L} + (SSHT) {swivel_step_half_turn_L} + (Sunanda) {sunanda_L} + (Balancing) {balancing_step_L} + (Figure 8) {figure_8_L} + (Barrel) {barrel_turn_L} = {total_L}")

# Print the equation for Back turns
print(f"Total Back Turns       = (Swivel Step) {swivel_step_B} + (SSHT) {swivel_step_half_turn_B} + (Sunanda) {sunanda_B} + (Balancing) {balancing_step_B} + (Figure 8) {figure_8_B} + (Barrel) {barrel_turn_B} = {total_B}")

print(f"\nFinal Answer: The dancer turns her/his right side {total_R} times, left side {total_L} times, and back {total_B} times.")