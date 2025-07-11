# My optimal first move is the center. The computer has an 8/8 chance of picking an open square.
# There are two symmetric cases for the computer's first move.

# Case 1: The computer plays on an edge square (4 of 8 possibilities).
p_comp_edge = "4/8"
# My analysis shows I can force a win if the computer makes an edge move.
p_win_if_edge = 1

# Case 2: The computer plays on a corner square (4 of 8 possibilities).
p_comp_corner = "4/8"
# The probability calculation for this case is more complex.
# P(Win|Corner) = P(Win on my 3rd move) + P(O blocks 1) * P(Win on my 4th move) + P(O blocks 1&2) * P(Win on my 5th move)
p_win_if_corner_calc = "(4/5) + (1/5)*(2/3) + (1/5)*(1/3)*(1/2)"
# Calculating the fraction
# 4/5 + 2/15 + 1/30 = 24/30 + 4/30 + 1/30 = 29/30
p_win_if_corner_val = "29/30"

# Total probability is the weighted average of the two cases.
final_calc_str = f"(({p_comp_corner}) * ({p_win_if_corner_val})) + (({p_comp_edge}) * {p_win_if_edge})"
final_result_num = (4/8) * (29/30) + (4/8) * 1
# (1/2)*(29/30) + (1/2)*1 = 29/60 + 30/60 = 59/60
final_result_frac = "59/60"

print("The maximum chance of winning is found by starting in the center.")
print("The computer can then make a move on a corner or an edge.")
print("\nCase 1: Computer plays on a corner (Probability = 4/8)")
print("My win probability in this case is calculated as the sum of win chances at each of the computer's subsequent moves:")
print(f"P(Win | Computer plays Corner) = {p_win_if_corner_calc} = {p_win_if_corner_val}")
print("\nCase 2: Computer plays on an edge (Probability = 4/8)")
print(f"My win probability in this case is 1 (a forced win).")
print("\nTotal Maximum Win Probability Calculation:")
print(f"P(Win) = (P(Corner) * P(Win|Corner)) + (P(Edge) * P(Win|Edge))")
print(f"P(Win) = {final_calc_str}")
print(f"P(Win) = (1/2) * (29/30) + (1/2) * 1 = 29/60 + 30/60 = {final_result_frac}")

print(f"\nThe maximum chance of winning is {final_result_frac}.")
