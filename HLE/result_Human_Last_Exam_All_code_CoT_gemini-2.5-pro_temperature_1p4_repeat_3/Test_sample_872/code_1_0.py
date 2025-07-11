import fractions

# Probabilities of the computer's first move after I take a corner.
# There are 8 possible moves for the computer.
p_comp_center = fractions.Fraction(1, 8)
p_comp_edge = fractions.Fraction(4, 8)
p_comp_adj_corner = fractions.Fraction(2, 8)
p_comp_opp_corner = fractions.Fraction(1, 8)

# My conditional win probabilities for each of the computer's moves.
# These are calculated based on my optimal second move and the computer's subsequent random move.
p_win_given_comp_center = fractions.Fraction(2, 6)      # I win if comp's 2nd move is a corner.
p_win_given_comp_edge = fractions.Fraction(1, 1)        # I have a forced win.
p_win_given_comp_adj_corner = fractions.Fraction(5, 6)  # I win if comp fails to make the single correct block.
p_win_given_comp_opp_corner = fractions.Fraction(4, 6)    # I win if comp's 2nd move is an edge.

# Calculate the contribution of each case to the total win probability.
win_prob_A = p_win_given_comp_center * p_comp_center
win_prob_B = p_win_given_comp_edge * p_comp_edge
win_prob_C = p_win_given_comp_adj_corner * p_comp_adj_corner
win_prob_D = p_win_given_comp_opp_corner * p_comp_opp_corner

# The total probability is the sum of the probabilities of these disjoint events.
total_win_prob = win_prob_A + win_prob_B + win_prob_C + win_prob_D

print("My optimal strategy is to start in a corner (or the center).")
print("The calculation for starting in a corner is as follows:\n")
print("P(Win) = P(Win|Comp plays Center) * P(Comp plays Center) +")
print("         P(Win|Comp plays Edge) * P(Comp plays Edge) +")
print("         P(Win|Comp plays Adj Corner) * P(Comp plays Adj Corner) +")
print("         P(Win|Comp plays Opp Corner) * P(Comp plays Opp Corner)\n")

print(f"P(Win) = ({p_win_given_comp_center.numerator}/{p_win_given_comp_center.denominator}) * ({p_comp_center.numerator}/{p_comp_center.denominator}) +")
print(f"         ({p_win_given_comp_edge.numerator}/{p_win_given_comp_edge.denominator}) * ({p_comp_edge.numerator}/{p_comp_edge.denominator}) +")
print(f"         ({p_win_given_comp_adj_corner.numerator}/{p_win_given_comp_adj_corner.denominator}) * ({p_comp_adj_corner.numerator}/{p_comp_adj_corner.denominator}) +")
print(f"         ({p_win_given_comp_opp_corner.numerator}/{p_win_given_comp_opp_corner.denominator}) * ({p_comp_opp_corner.numerator}/{p_comp_opp_corner.denominator})\n")

print(f"P(Win) = {win_prob_A} + {win_prob_B} + {win_prob_C} + {win_prob_D}\n")
print(f"P(Win) = {win_prob_A.limit_denominator()} + {win_prob_B.limit_denominator()} + {win_prob_C.limit_denominator()} + {win_prob_D.limit_denominator()}\n")
# To show the sum with a common denominator
common_denominator = 48
print(f"P(Win) = {win_prob_A.numerator * (common_denominator // win_prob_A.denominator)}/{common_denominator} + {win_prob_B.numerator * (common_denominator // win_prob_B.denominator)}/{common_denominator} + {win_prob_C.numerator * (common_denominator // win_prob_C.denominator)}/{common_denominator} + {win_prob_D.numerator * (common_denominator // win_prob_D.denominator)}/{common_denominator}\n")
print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")