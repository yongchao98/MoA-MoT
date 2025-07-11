from fractions import Fraction

# Step 1: Explain the optimal strategy.
print("To find the maximum chance of winning, we must use an optimal strategy.")
print("The optimal first move in Tic Tac Toe is to place your 'X' in a corner.")
print("After we play in a corner, the computer, playing randomly, places 'O' in one of the 8 remaining squares.")
print("\nWe analyze the computer's move based on whether it plays in the center or not.")

# Step 2: Calculate the win probability when the computer does NOT play center.
prob_o_not_center = Fraction(7, 8)
win_prob_if_not_center = Fraction(1, 1)
contribution_not_center = prob_o_not_center * win_prob_if_not_center

print(f"\nCase 1: The computer does NOT play in the center.")
print(f"The probability of this is 7/8.")
print("In this scenario, our optimal second move guarantees a win against a random opponent.")
print(f"Thus, our win probability is 1.")
print(f"Contribution to total win P: ({prob_o_not_center.numerator}/{prob_o_not_center.denominator}) * {win_prob_if_not_center.numerator}/{win_prob_if_not_center.denominator} = {contribution_not_center.numerator}/{contribution_not_center.denominator}")


# Step 3: Calculate the win probability when the computer plays center.
prob_o_center = Fraction(1, 8)
print(f"\nCase 2: The computer plays in the center.")
print(f"The probability of this is 1/8.")
print("Our best response is the opposite corner. The computer then has 6 squares for its second move.")
# Win probability in this sub-case is derived from the computer's 2nd move.
# P(win | O in center) = P(O2 is edge) * P(win|O2 is edge) + P(O2 is corner) * P(win|O2 is corner)
# = (4/6 * 1) + (2/6 * 3/4) = 2/3 + 1/4 = 11/12
win_prob_if_center = Fraction(11, 12)
contribution_center = prob_o_center * win_prob_if_center
print(f"In this scenario, a detailed breakdown shows our win probability is 11/12.")
print(f"Contribution to total win P: ({prob_o_center.numerator}/{prob_o_center.denominator}) * ({win_prob_if_center.numerator}/{win_prob_if_center.denominator}) = {contribution_center.numerator}/{contribution_center.denominator}")

# Step 4: Calculate the total win probability.
total_win_prob = contribution_not_center + contribution_center
c_not_center_num = contribution_not_center.numerator
c_not_center_den = contribution_not_center.denominator
c_center_num = contribution_center.numerator
c_center_den = contribution_center.denominator
total_num = total_win_prob.numerator
total_den = total_win_prob.denominator

# To show the common denominator calculation
c_not_center_num_adj = c_not_center_num * (total_den // c_not_center_den)
c_center_num_adj = c_center_num * (total_den // c_center_den)

print("\nFinal Calculation:")
print(f"The total probability is the sum of the contributions from each case.")
print(f"Total P(Win) = {c_not_center_num}/{c_not_center_den} + {c_center_num}/{c_center_den}")
print(f"Using a common denominator: {c_not_center_num_adj}/{total_den} + {c_center_num_adj}/{total_den} = {total_num}/{total_den}")
