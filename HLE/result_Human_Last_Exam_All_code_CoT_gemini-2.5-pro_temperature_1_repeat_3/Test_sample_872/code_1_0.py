import fractions

# Step 1: Define the initial probabilities for the computer's first move.
# I play first in the center. The computer has 8 squares left.
# 4 are corners, 4 are edges.
prob_O_plays_corner = fractions.Fraction(4, 8)
prob_O_plays_edge = fractions.Fraction(4, 8)

print("My Optimal Strategy:")
print("1. I place my first 'X' in the center square.")
print("2. The computer randomly places an 'O' in one of the 8 remaining squares.")
print(f"   - The probability of the computer choosing a corner is {prob_O_plays_corner.numerator}/{prob_O_plays_corner.denominator}.")
print(f"   - The probability of the computer choosing an edge is {prob_O_plays_edge.numerator}/{prob_O_plays_edge.denominator}.")
print("-" * 20)

# Step 2: Analyze Case 1: Computer plays a corner.
# My best response is to place my second 'X' on an adjacent edge, creating a threat.
# The computer has 6 squares left for its second move.
prob_O2_blocks_threat1 = fractions.Fraction(1, 6)
prob_O2_not_block_threat1 = fractions.Fraction(5, 6)

# If O2 blocks, I place X3 to create a new threat. The computer (O3) has 4 squares left.
# My move creates a situation where the computer must block a specific square to prevent a loss.
# Of the 4 available moves for the computer, 3 lead to me winning, and 1 leads to a draw (which is a loss for me).
win_prob_if_O2_blocks = fractions.Fraction(3, 4)

# The win probability if the computer's first move was a corner:
win_prob_given_O_corner = (prob_O2_not_block_threat1 * 1) + (prob_O2_blocks_threat1 * win_prob_if_O2_blocks)

print("Analysis for Case 1 (Computer plays a corner):")
print("  - I place my second 'X' on an adjacent edge, creating a threat.")
print(f"  - The computer fails to block my threat with probability {prob_O2_not_block_threat1.numerator}/{prob_O2_not_block_threat1.denominator}, and I win.")
print(f"  - The computer blocks my threat with probability {prob_O2_blocks_threat1.numerator}/{prob_O2_blocks_threat1.denominator}.")
print(f"    - If it blocks, my next move forces a situation where I win with probability {win_prob_if_O2_blocks.numerator}/{win_prob_if_O2_blocks.denominator}.")
print(f"  - Total win probability in this case: ({prob_O2_not_block_threat1.numerator}/{prob_O2_not_block_threat1.denominator}) * 1 + ({prob_O2_blocks_threat1.numerator}/{prob_O2_blocks_threat1.denominator}) * ({win_prob_if_O2_blocks.numerator}/{win_prob_if_O2_blocks.denominator}) = {win_prob_given_O_corner.numerator}/{win_prob_given_O_corner.denominator}")
print("-" * 20)

# Step 3: Analyze Case 2: Computer plays an edge.
# My best response is to place my second 'X' in an adjacent corner, creating a threat.
# The computer has 6 squares left for its second move.
prob_O2_blocks_threat2 = fractions.Fraction(1, 6)
prob_O2_not_block_threat2 = fractions.Fraction(5, 6)

# If the computer blocks, my next move creates a "fork" (two ways to win), which is an unstoppable attack.
win_prob_if_O2_blocks_edge = 1

# The win probability if the computer's first move was an edge:
win_prob_given_O_edge = (prob_O2_not_block_threat2 * 1) + (prob_O2_blocks_threat2 * win_prob_if_O2_blocks_edge)

print("Analysis for Case 2 (Computer plays an edge):")
print("  - I place my second 'X' on an adjacent corner, creating a threat.")
print(f"  - The computer fails to block with probability {prob_O2_not_block_threat2.numerator}/{prob_O2_not_block_threat2.denominator}, and I win.")
print(f"  - The computer blocks with probability {prob_O2_blocks_threat2.numerator}/{prob_O2_blocks_threat2.denominator}.")
print(f"    - If it blocks, my next move creates a fork, guaranteeing a win (probability 1/1).")
print(f"  - Total win probability in this case: ({prob_O2_not_block_threat2.numerator}/{prob_O2_not_block_threat2.denominator}) * 1 + ({prob_O2_blocks_threat2.numerator}/{prob_O2_blocks_threat2.denominator}) * 1 = {win_prob_given_O_edge.numerator}/{win_prob_given_O_edge.denominator}")
print("-" * 20)

# Step 4: Calculate the total win probability.
total_win_prob = (prob_O_plays_corner * win_prob_given_O_corner) + (prob_O_plays_edge * win_prob_given_O_edge)

print("Final Calculation:")
print("The maximum chance of winning is the sum of probabilities of each case:")
print(f"P(Win) = P(O plays corner) * P(Win | O plays corner) + P(O plays edge) * P(Win | O plays edge)")
print(f"P(Win) = ({prob_O_plays_corner.numerator}/{prob_O_plays_corner.denominator}) * ({win_prob_given_O_corner.numerator}/{win_prob_given_O_corner.denominator}) + ({prob_O_plays_edge.numerator}/{prob_O_plays_edge.denominator}) * ({win_prob_given_O_edge.numerator}/{win_prob_given_O_edge.denominator})")
print(f"P(Win) = {prob_O_plays_corner * win_prob_given_O_corner} + {prob_O_plays_edge * win_prob_given_O_edge}")
print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")
print("-" * 20)
print(f"The maximum chance of winning is {total_win_prob.numerator}/{total_win_prob.denominator}.")
