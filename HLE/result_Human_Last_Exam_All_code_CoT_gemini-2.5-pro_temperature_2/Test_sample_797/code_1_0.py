# Define the number of times each move is performed.
# Even if a move has a count number (like 8 or 16), it is still performed once in the sequence unless specified otherwise.
swivel_step_reps = 1
swivel_step_half_turn_reps = 1
sunanda_reps = 1
balancing_step_reps = 1
figure_8_reps = 8 # This move itself does not involve a turn.
barrel_turn_reps = 1

# Define the number of right, left, and back turns for each type of move based on our logic.
# [Right, Left, Back]
swivel_step_turns = [0, 0, 0]
swivel_step_half_turn_turns = [0, 0, 1]
sunanda_turns = [1, 1, 0]
balancing_step_turns = [0, 0, 0]
figure_8_turn_per_rep = [0, 0, 0] # A standard figure 8 has no turn.
barrel_turn_turns = [1, 1, 0]

# Calculate the total number of turns for each side.
total_right_turns = (swivel_step_reps * swivel_step_turns[0] +
                     swivel_step_half_turn_reps * swivel_step_half_turn_turns[0] +
                     sunanda_reps * sunanda_turns[0] +
                     balancing_step_reps * balancing_step_turns[0] +
                     figure_8_reps * figure_8_turn_per_rep[0] +
                     barrel_turn_reps * barrel_turn_turns[0])

total_left_turns = (swivel_step_reps * swivel_step_turns[1] +
                    swivel_step_half_turn_reps * swivel_step_half_turn_turns[1] +
                    sunanda_reps * sunanda_turns[1] +
                    balancing_step_reps * balancing_step_turns[1] +
                    figure_8_reps * figure_8_turn_per_rep[1] +
                    barrel_turn_reps * barrel_turn_turns[1])

total_back_turns = (swivel_step_reps * swivel_step_turns[2] +
                    swivel_step_half_turn_reps * swivel_step_half_turn_turns[2] +
                    sunanda_reps * sunanda_turns[2] +
                    balancing_step_reps * balancing_step_turns[2] +
                    figure_8_reps * figure_8_turn_per_rep[2] +
                    barrel_turn_reps * barrel_turn_turns[2])

# Print the breakdown of the calculation for clarity.
print("Calculating total right side turns:")
print(f"Sunanda ({sunanda_reps}x) adds {sunanda_turns[0]} right turn(s).")
print(f"Barrel Turn ({barrel_turn_reps}x) adds {barrel_turn_turns[0]} right turn(s).")
print(f"Total Right Turns = {sunanda_turns[0]} + {barrel_turn_turns[0]} = {total_right_turns}\n")

print("Calculating total left side turns:")
print(f"Sunanda ({sunanda_reps}x) adds {sunanda_turns[1]} left turn(s).")
print(f"Barrel Turn ({barrel_turn_reps}x) adds {barrel_turn_turns[1]} left turn(s).")
print(f"Total Left Turns = {sunanda_turns[1]} + {barrel_turn_turns[1]} = {total_left_turns}\n")

print("Calculating total back turns:")
print(f"Swivel Step Half Turn ({swivel_step_half_turn_reps}x) adds {swivel_step_half_turn_turns[2]} back turn(s).")
print(f"Total Back Turns = {swivel_step_half_turn_turns[2]} = {total_back_turns}\n")

print(f"Final Answer: The dancer turns her/his right side {total_right_turns} times, left side {total_left_turns} times, and back {total_back_turns} time(s).")
print(f"This corresponds to the ratio {total_right_turns}/{total_left_turns}/{total_back_turns}.")