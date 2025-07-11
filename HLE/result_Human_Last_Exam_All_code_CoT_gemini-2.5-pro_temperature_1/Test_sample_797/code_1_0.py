# Initialize counts for each direction
right_turns = 0
left_turns = 0
back_turns = 0

# --- Breakdown of the dance sequence ---

# 1. Swivel Step (8 counts)
# Interpreted as a quarter turn right and a quarter turn left.
swivel_step_r = 1
swivel_step_l = 1
swivel_step_b = 0
right_turns += swivel_step_r
left_turns += swivel_step_l
back_turns += swivel_step_b

# 2. Swivel Step Half Turn (16 counts)
# A 180-degree turn and return shows both sides.
swivel_half_turn_r = 1
swivel_half_turn_l = 1
swivel_half_turn_b = 0 # Back turn is counted in the more deliberate Sunanda
right_turns += swivel_half_turn_r
left_turns += swivel_half_turn_l
back_turns += swivel_half_turn_b

# 3. Sunanda (once)
# Two 180-degree turns, this is the move where the back turn is counted.
sunanda_r = 1
sunanda_l = 1
sunanda_b = 1
right_turns += sunanda_r
left_turns += sunanda_l
back_turns += sunanda_b

# 4. Balancing Step (once)
# No turns.
balancing_step_r = 0
balancing_step_l = 0
balancing_step_b = 0
right_turns += balancing_step_r
left_turns += balancing_step_l
back_turns += balancing_step_b

# 5. Figure 8 (8 times)
# No turns.
figure_8_r = 0
figure_8_l = 0
figure_8_b = 0
right_turns += figure_8_r
left_turns += figure_8_l
back_turns += figure_8_b

# 6. Barrel Turn (once)
# A full 360-degree spin shows both sides.
barrel_turn_r = 1
barrel_turn_l = 1
barrel_turn_b = 0 # Back turn is counted in Sunanda
right_turns += barrel_turn_r
left_turns += barrel_turn_l
back_turns += barrel_turn_b

# --- Final Calculation ---
print("Calculating the total number of turns for right side, left side, and back.")

print("\nRight side turns:")
print(f"Swivel Step ({swivel_step_r}) + Swivel Step Half Turn ({swivel_half_turn_r}) + Sunanda ({sunanda_r}) + Balancing Step ({balancing_step_r}) + Figure 8 ({figure_8_r}) + Barrel Turn ({barrel_turn_r}) = {right_turns}")

print("\nLeft side turns:")
print(f"Swivel Step ({swivel_step_l}) + Swivel Step Half Turn ({swivel_half_turn_l}) + Sunanda ({sunanda_l}) + Balancing Step ({balancing_step_l}) + Figure 8 ({figure_8_l}) + Barrel Turn ({barrel_turn_l}) = {left_turns}")

print("\nBack turns:")
print(f"Swivel Step ({swivel_step_b}) + Swivel Step Half Turn ({swivel_half_turn_b}) + Sunanda ({sunanda_b}) + Balancing Step ({balancing_step_b}) + Figure 8 ({figure_8_b}) + Barrel Turn ({barrel_turn_b}) = {back_turns}")

print(f"\nTotal turns are: Right/{right_turns}, Left/{left_turns}, Back/{back_turns}")