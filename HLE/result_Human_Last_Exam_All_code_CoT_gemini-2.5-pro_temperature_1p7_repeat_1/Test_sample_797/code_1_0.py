# Plan:
# 1. Define variables representing the number of right, left, and back turns for each move in the sequence.
# 2. Calculate the total for each orientation by summing the turns from each performed move.
# 3. Print the final calculation for each orientation, showing how the total was derived from each move.

# Turns from each move (Right, Left, Back)
swivel_step = (1, 0, 0)
swivel_step_half_turn = (1, 0, 1)
sunanda = (1, 1, 0)
balancing_step = (0, 0, 0)
figure_8 = (0, 0, 0) # This move is done 8 times, but has 0 turns.
barrel_turn = (1, 1, 0)

# Calculate the total for Right side turns
total_right = swivel_step[0] + swivel_step_half_turn[0] + sunanda[0] + balancing_step[0] + figure_8[0] + barrel_turn[0]
print(f"Right Side Turns = {swivel_step[0]} (Swivel Step) + {swivel_step_half_turn[0]} (Swivel Step Half Turn) + {sunanda[0]} (Sunanda) + {balancing_step[0]} (Balancing Step) + {figure_8[0]} (Figure 8) + {barrel_turn[0]} (Barrel Turn) = {total_right}")

# Calculate the total for Left side turns
total_left = swivel_step[1] + swivel_step_half_turn[1] + sunanda[1] + balancing_step[1] + figure_8[1] + barrel_turn[1]
print(f"Left Side Turns = {swivel_step[1]} (Swivel Step) + {swivel_step_half_turn[1]} (Swivel Step Half Turn) + {sunanda[1]} (Sunanda) + {balancing_step[1]} (Balancing Step) + {figure_8[1]} (Figure 8) + {barrel_turn[1]} (Barrel Turn) = {total_left}")

# Calculate the total for Back turns
total_back = swivel_step[2] + swivel_step_half_turn[2] + sunanda[2] + balancing_step[2] + figure_8[2] + barrel_turn[2]
print(f"Back Turns = {swivel_step[2]} (Swivel Step) + {swivel_step_half_turn[2]} (Swivel Step Half Turn) + {sunanda[2]} (Sunanda) + {balancing_step[2]} (Balancing Step) + {figure_8[2]} (Figure 8) + {barrel_turn[2]} (Barrel Turn) = {total_back}")

print(f"\nFinal count: Right/{total_right}, Left/{total_left}, Back/{total_back}")