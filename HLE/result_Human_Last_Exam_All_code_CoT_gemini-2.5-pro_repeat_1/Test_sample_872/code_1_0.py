from fractions import Fraction

# Step 1: Define the probabilities of the computer's first move based on my optimal start in the center.
# There are 4 corner squares and 4 edge squares out of 8 total available squares.
prob_O1_corner = Fraction(4, 8)
prob_O1_edge = Fraction(4, 8)

# Step 2: Calculate the win probability for Case A (computer plays a corner).
# My best response (X2) is to take another corner, creating a threat.
# The computer has 6 squares left for its move (O2) and must block one specific square.
prob_O2_fails_to_block_threat1 = Fraction(5, 6)
prob_O2_blocks_threat1 = Fraction(1, 6)

# If the computer blocks (the 1/6 case), I make another move (X3) to create a new threat.
# Now there are 4 squares left for the computer's move (O3).
prob_O3_fails_to_block_threat2 = Fraction(3, 4)

# My win probability in the "corner" branch is the sum of probabilities of me winning at different stages.
# I win if O2 fails to block, OR if O2 blocks AND O3 fails to block.
win_prob_if_corner = prob_O2_fails_to_block_threat1 + prob_O2_blocks_threat1 * prob_O3_fails_to_block_threat2

# Step 3: Calculate the win probability for Case B (computer plays an edge).
# My best response (X2) creates a threat.
# The computer has 6 squares left for O2.
prob_O2_fails_to_block_fork_setup = Fraction(5, 6)
prob_O2_blocks_fork_setup = Fraction(1, 6)

# If the computer blocks (1/6 chance), my next move (X3) forces a block from the computer
# and simultaneously creates a "fork" (two ways to win). The computer can only block
# one path, so my win is guaranteed.
win_prob_if_block_leads_to_fork = Fraction(1, 1)

# My win probability in the "edge" branch. I win if O2 fails to block, OR if O2 blocks (which leads to my guaranteed win).
win_prob_if_edge = prob_O2_fails_to_block_fork_setup * win_prob_if_block_leads_to_fork + prob_O2_blocks_fork_setup * win_prob_if_block_leads_to_fork


# Step 4: Calculate the total maximum win probability using the law of total probability.
total_win_prob = win_prob_if_corner * prob_O1_corner + win_prob_if_edge * prob_O1_edge

# Step 5: Print out the steps of the calculation.
print("Here is the step-by-step calculation for the maximum chance of winning:")
print("\nMy optimal first move is to take the center square.")
print("The computer then plays a random unfilled square. There are two cases for its first move:")
print(f"1. The computer plays a corner ({prob_O1_corner.numerator}/{prob_O1_corner.denominator} probability).")
print(f"2. The computer plays an edge ({prob_O1_edge.numerator}/{prob_O1_edge.denominator} probability).")

print("\n--- Analysis of Case 1 (Computer plays a corner) ---")
print("My best response creates a threat. The computer must block it on its next move.")
print(f"The computer has 6 open squares. The probability it fails to block is {prob_O2_fails_to_block_threat1.numerator}/{prob_O2_fails_to_block_threat1.denominator}.")
print(f"If it fails, I win. If it successfully blocks (a {prob_O2_blocks_threat1.numerator}/{prob_O2_blocks_threat1.denominator} chance), I create a new threat.")
print(f"The computer now has 4 open squares. The probability it fails to block this second threat is {prob_O3_fails_to_block_threat2.numerator}/{prob_O3_fails_to_block_threat2.denominator}.")
print("The total win probability in this case is: P(Win|Corner) = P(O fails 1st time) + P(O blocks 1st time) * P(O fails 2nd time)")
print(f"P(Win|Corner) = {prob_O2_fails_to_block_threat1.numerator}/{prob_O2_fails_to_block_threat1.denominator} + {prob_O2_blocks_threat1.numerator}/{prob_O2_blocks_threat1.denominator} * {prob_O3_fails_to_block_threat2.numerator}/{prob_O3_fails_to_block_threat2.denominator} = {win_prob_if_corner.numerator}/{win_prob_if_corner.denominator}")

print("\n--- Analysis of Case 2 (Computer plays an edge) ---")
print("My best response creates a threat that, if blocked, leads to a guaranteed win for me (a 'fork').")
print(f"The computer has 6 open squares. If it fails to block (a {prob_O2_fails_to_block_fork_setup.numerator}/{prob_O2_fails_to_block_fork_setup.denominator} chance), I win.")
print(f"If it blocks (a {prob_O2_blocks_fork_setup.numerator}/{prob_O2_blocks_fork_setup.denominator} chance), my next move creates a guaranteed fork.")
print(f"A guaranteed win means the probability of winning is {win_prob_if_block_leads_to_fork.numerator}/{win_prob_if_block_leads_to_fork.denominator}.")
print("The total win probability in this case is: P(Win|Edge) = P(O fails to block) * 1 + P(O blocks) * 1")
print(f"P(Win|Edge) = {prob_O2_fails_to_block_fork_setup.numerator}/{prob_O2_fails_to_block_fork_setup.denominator} + {prob_O2_blocks_fork_setup.numerator}/{prob_O2_blocks_fork_setup.denominator} = {win_prob_if_edge.numerator}/{win_prob_if_edge.denominator}")

print("\n--- Final Calculation ---")
print("Total P(Win) = P(Win|Corner) * P(Corner) + P(Win|Edge) * P(Edge)")
print(f"Total P(Win) = ({win_prob_if_corner.numerator}/{win_prob_if_corner.denominator}) * ({prob_O1_corner.numerator}/{prob_O1_corner.denominator}) + ({win_prob_if_edge.numerator}/{win_prob_if_edge.denominator}) * ({prob_O1_edge.numerator}/{prob_O1_edge.denominator})")
final_term1 = win_prob_if_corner * prob_O1_corner
final_term2 = win_prob_if_edge * prob_O1_edge
print(f"Total P(Win) = {final_term1.numerator}/{final_term1.denominator} + {final_term2.numerator}/{final_term2.denominator}")
print(f"Total P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")