from fractions import Fraction

# The problem is a probability calculation based on an optimal strategy.
# My optimal first move is to take the center square.
# The computer's first move is random, with 8 possibilities.
# We can group the computer's move into two symmetric cases.

# Case 1: The computer plays an edge square (4 of 8 squares).
prob_comp_edge = Fraction(4, 8)
# If the computer plays an edge, my optimal response guarantees a win.
# For example: Me(1,1), Comp(0,1), Me(0,0). Comp must not play (2,2) to block my
# immediate threat. If it does, I create a fork and win. If it doesn't, I win.
# Therefore, my win probability in this case is 1.
prob_win_if_comp_edge = Fraction(1, 1)
contribution_from_edge_case = prob_comp_edge * prob_win_if_comp_edge

# Case 2: The computer plays a corner square (4 of 8 squares).
prob_comp_corner = Fraction(4, 8)
# If the computer plays a corner, my optimal response is the opposite corner.
# e.g., Me(1,1), Comp(0,0), Me(2,2).
# Now, the computer has 6 squares left for its second move.
# We analyze the sub-cases for the computer's second move.

#   Sub-case 2a: Computer's second move creates a threat (e.g., at (0,1) or (1,0)).
#   This happens with probability 2/6.
prob_comp_threatens_in_corner_case = Fraction(2, 6)
#   I must block. After my block, I have a threat the computer must block.
#   The computer has 4 moves, and only 1 is correct to force a tie.
#   So, the computer fails to block with probability 3/4, and I win.
prob_win_if_comp_threatens = Fraction(3, 4)

#   Sub-case 2b: Computer's second move does not create a threat.
#   This happens with probability 4/6.
prob_comp_no_threat_in_corner_case = Fraction(4, 6)
#   In this case, I can make a move that creates a fork, guaranteeing a win.
#   So, my win probability is 1.
prob_win_if_comp_no_threat = Fraction(1, 1)

# The total win probability for Case 2 is the weighted average of the sub-cases.
prob_win_if_comp_corner = (prob_comp_threatens_in_corner_case * prob_win_if_comp_threatens) + \
                          (prob_comp_no_threat_in_corner_case * prob_win_if_comp_no_threat)
contribution_from_corner_case = prob_comp_corner * prob_win_if_comp_corner

# The total maximum win probability is the sum of the contributions from both cases.
total_win_prob = contribution_from_edge_case + contribution_from_corner_case

# --- Output the explanation and calculation ---
print("To maximize the chance of winning, the optimal first move is to take the center square.")
print("We then analyze the computer's random response, which falls into two main cases:")

print("\nCase 1: Computer plays an edge square.")
print(f"  - Probability of this case: {prob_comp_edge.numerator}/{prob_comp_edge.denominator}")
print(f"  - Your win probability with optimal play: {prob_win_if_comp_edge.numerator}/{prob_win_if_comp_edge.denominator}")

print("\nCase 2: Computer plays a corner square.")
print(f"  - Probability of this case: {prob_comp_corner.numerator}/{prob_comp_corner.denominator}")
print("  - Your win probability here is calculated from its sub-cases:")
print(f"    P(Win|Corner) = P(Comp Threatens) * P(Win|Threat) + P(Comp No Threat) * P(Win|No Threat)")
print(f"    P(Win|Corner) = ({prob_comp_threatens_in_corner_case}) * ({prob_win_if_comp_threatens}) + ({prob_comp_no_threat_in_corner_case}) * ({prob_win_if_comp_no_threat})")
print(f"    P(Win|Corner) = {prob_win_if_comp_corner.numerator}/{prob_win_if_comp_corner.denominator}")

print("\nCombining these cases gives the total maximum win probability.")
print("The final equation is:")
print(f"P(Total Win) = P(Comp Edge) * P(Win|Edge) + P(Comp Corner) * P(Win|Corner)")
print(f"P(Total Win) = ({prob_comp_edge}) * ({prob_win_if_comp_edge}) + ({prob_comp_corner}) * ({prob_win_if_comp_corner})")
print(f"P(Total Win) = ({contribution_from_edge_case}) + ({contribution_from_corner_case})")
print(f"P(Total Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")

print(f"\nThus, the maximum chance of winning you can give yourself is {total_win_prob.numerator}/{total_win_prob.denominator}.")

<<<23/24>>>