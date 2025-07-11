import fractions

# As derived from the random walk model, the probability of Theo winning a single game is 1/6.
prob_theo_wins = fractions.Fraction(1, 6)

# The event "Theo wins for the first time only after at least five games"
# is equivalent to "Theo does not win in any of the first four games".
# The probability of Theo not winning a single game is 1 - P(Theo wins).
prob_theo_not_wins = 1 - prob_theo_wins

# The number of games in which Theo does not win.
num_games_no_win = 4

# The final probability is P(Theo does not win) raised to the power of 4,
# as the games are independent events.
final_prob = prob_theo_not_wins ** num_games_no_win

# Extract numbers for the equation string
numerator_base = prob_theo_not_wins.numerator
denominator_base = prob_theo_not_wins.denominator
power = num_games_no_win
final_numerator = final_prob.numerator
final_denominator = final_prob.denominator

# Print the final result and the equation used to calculate it.
print("The probability of Theo winning a single game is 1/6.")
print("Therefore, the probability of Theo NOT winning a single game is 1 - 1/6 = 5/6.")
print("\nThe probability that Theo wins for the first time only after at least five games is the probability that he does not win the first four games.")
print(f"\nThis is calculated as: ({numerator_base}/{denominator_base})^{power}")
print(f"Final Calculation: {numerator_base}^{power} / {denominator_base}^{power} = {final_numerator} / {final_denominator}")
print(f"\nThe final probability as a fraction is: {final_numerator}/{final_denominator}")
print(f"The final probability as a decimal is approximately: {float(final_prob):.6f}")