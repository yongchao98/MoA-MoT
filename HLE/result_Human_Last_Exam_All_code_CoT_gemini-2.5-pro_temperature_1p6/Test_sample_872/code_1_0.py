from fractions import Fraction

# Define the initial probabilities for the computer's first move
prob_o_plays_corner = Fraction(4, 8)
prob_o_plays_edge = Fraction(4, 8)

# --- Case 1: Computer plays a corner ---
# Probability computer fails to make the first block
prob_win_on_turn3_corner = Fraction(5, 6)

# Probability computer makes the first block
prob_o_blocks_turn2_corner = Fraction(1, 6)
# If the computer blocks, a deeper analysis shows the subsequent win probability for the player is 7/8.
# This is derived from: 3/4 + (1/4 * 1/2) = 7/8
prob_win_after_block_corner = Fraction(7, 8)

# Total win probability if computer chose a corner
prob_win_if_corner = prob_win_on_turn3_corner + (prob_o_blocks_turn2_corner * prob_win_after_block_corner)

# --- Case 2: Computer plays an edge ---
# Probability computer fails to make the first block
prob_win_on_turn3_edge = Fraction(5, 6)

# Probability computer makes the first block
prob_o_blocks_turn2_edge = Fraction(1, 6)
# If the computer blocks, the player can force a fork, guaranteeing a win.
prob_win_after_block_edge = Fraction(1, 1)

# Total win probability if computer chose an edge
prob_win_if_edge = prob_win_on_turn3_edge + (prob_o_blocks_turn2_edge * prob_win_after_block_edge)

# --- Final Calculation ---
# The total probability is the weighted average of the two cases.
total_prob = (prob_o_plays_corner * prob_win_if_corner) + (prob_o_plays_edge * prob_win_if_edge)

print("This script calculates the maximum probability of winning Tic Tac Toe against a random opponent.")
print("\nStep 1: Calculate the win probability if the computer plays a corner first.")
print(f"P(Win | O plays corner) = P(O doesn't block) + P(O blocks) * P(Win after block)")
print(f"P(Win | O plays corner) = {prob_win_on_turn3_corner} + {prob_o_blocks_turn2_corner} * {prob_win_after_block_corner} = {prob_win_if_corner}")

print("\nStep 2: Calculate the win probability if the computer plays an edge first.")
print(f"P(Win | O plays edge) = P(O doesn't block) + P(O blocks) * P(Guaranteed Win)")
print(f"P(Win | O plays edge) = {prob_win_on_turn3_edge} + {prob_o_blocks_turn2_edge} * {prob_win_after_block_edge} = {prob_win_if_edge}")

print("\nStep 3: Combine the probabilities for the final answer.")
print(f"P(Total) = P(O plays corner) * P(Win | O plays corner) + P(O plays edge) * P(Win | O plays edge)")
print(f"P(Total) = {prob_o_plays_corner} * {prob_win_if_corner} + {prob_o_plays_edge} * {prob_win_if_edge}")
print(f"P(Total) = {prob_o_plays_corner * prob_win_if_corner} + {prob_o_plays_edge * prob_win_if_edge}")
print(f"P(Total) = {total_prob}")

# Final Answer
print("\nThe maximum chance of winning is:")
print(f"<<<{total_prob}>>>")
