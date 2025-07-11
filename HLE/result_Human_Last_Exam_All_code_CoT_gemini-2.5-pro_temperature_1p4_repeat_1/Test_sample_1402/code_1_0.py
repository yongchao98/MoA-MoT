import fractions

# Step 1: Calculate the probability of Theo winning a single game (P_TW).
# This is based on a random walk model where D = H - T.
# The game ends if D hits 0 (draw), 3 (Alex wins), or -3 (Theo wins).
# Let p_T(d) be the probability Theo wins starting from state d.
# We have boundary conditions: p_T(-3)=1, p_T(3)=0, p_T(0)=0.
# The recurrence is p_T(d) = 0.5 * p_T(d+1) + 0.5 * p_T(d-1).

# Solving for p_T(-1) and p_T(1):
# For d > 0:
# p_T(1) = 0.5*p_T(2) + 0.5*p_T(0) = 0.5*p_T(2)
# p_T(2) = 0.5*p_T(3) + 0.5*p_T(1) = 0.5*0 + 0.5*p_T(1) = 0.5*p_T(1)
# Substituting gives p_T(1) = 0.5 * (0.5 * p_T(1)) = 0.25 * p_T(1), which means p_T(1) = 0.
p_T_from_state_1 = fractions.Fraction(0)

# For d < 0:
# p_T(-1) = 0.5*p_T(0) + 0.5*p_T(-2) = 0.5*p_T(-2)
# p_T(-2) = 0.5*p_T(-1) + 0.5*p_T(-3) = 0.5*p_T(-1) + 0.5*1
# Substituting gives p_T(-2) = 0.5 * (0.5 * p_T(-2)) + 0.5, so 0.75*p_T(-2) = 0.5.
# This yields p_T(-2) = 2/3, and therefore p_T(-1) = 0.5 * (2/3) = 1/3.
p_T_from_state_neg_1 = fractions.Fraction(1, 3)

# The game starts with one toss (H or T), so we average the probabilities from states 1 and -1.
prob_theo_wins = fractions.Fraction(1, 2) * p_T_from_state_1 + fractions.Fraction(1, 2) * p_T_from_state_neg_1

# Step 2: Calculate the probability of Theo NOT winning a single game.
prob_theo_not_win = 1 - prob_theo_wins

# Step 3: Calculate the probability that Theo does not win in the first 4 games.
# This is equivalent to "Theo wins for the first time only after at least five games".
num_games_without_win = 4
final_prob = prob_theo_not_win ** num_games_without_win

# Step 4: Output the final equation and the result.
numerator = prob_theo_not_win.numerator
denominator = prob_theo_not_win.denominator
exponent = num_games_without_win
result_numerator = final_prob.numerator
result_denominator = final_prob.denominator

print(f"The probability of Theo winning a single game is {prob_theo_wins}.")
print(f"The probability of Theo not winning a single game is {prob_theo_not_win}.")
print("\nThe probability that Theo wins for the first time only after at least five games is the probability that he does not win in the first four games.")
print("\nThe final calculation is:")
print(f"({numerator}/{denominator})^({exponent}) = {result_numerator}/{result_denominator}")