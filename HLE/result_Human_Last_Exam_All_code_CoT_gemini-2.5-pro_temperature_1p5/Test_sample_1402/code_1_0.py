import fractions

# Step 1: Define the probabilities for a single game.
# Based on a random walk analysis, we determined the probabilities of each outcome.
p_theo_win_num = 1
p_theo_win_den = 6
p_theo_win = fractions.Fraction(p_theo_win_num, p_theo_win_den)

print("This problem involves a sequence of independent games. Let's first analyze a single game.")
print(f"The probability of Theo winning a single game (P_T) is {p_theo_win.numerator}/{p_theo_win.denominator}.")
print(f"The probability of Alex winning a single game (P_A) is also {p_theo_win.numerator}/{p_theo_win.denominator} by symmetry.")
print(f"The probability of a draw (P_D) is 1 - P_T - P_A = 1 - 1/6 - 1/6 = 2/3.")
print("-" * 50)

# Step 2: Interpret the question about the sequence of games.
print("The question asks for the probability that 'Theo wins for the first time only after at least five games'.")
print("This means Theo does not win in game 1, game 2, game 3, or game 4.")
print("-" * 50)

# Step 3: Calculate the probability that Theo does not win a single game.
p_not_theo_win = 1 - p_theo_win
print("The probability of Theo NOT winning a single game is 1 - P_T.")
print(f"P(Not T) = 1 - {p_theo_win.numerator}/{p_theo_win.denominator} = {p_not_theo_win.numerator}/{p_not_theo_win.denominator}")
print("-" * 50)

# Step 4: Calculate the final probability for the first four games.
num_games = 4
final_prob = p_not_theo_win ** num_games

print(f"Since the games are independent, we calculate the probability that Theo does not win in the first {num_games} games.")
print("The calculation is P(Not T) * P(Not T) * P(Not T) * P(Not T) = P(Not T)^4.")
print("\nFinal Equation:")

# Build and print the equation string with each term
equation_terms = [f"({p_not_theo_win.numerator}/{p_not_theo_win.denominator})" for _ in range(num_games)]
print(f"{' * '.join(equation_terms)} = ({p_not_theo_win.numerator}/{p_not_theo_win.denominator})^{num_games}")

# Print the calculation of the numerator and denominator
print(f"= {p_not_theo_win.numerator}^{num_games} / {p_not_theo_win.denominator}^{num_games}")
print(f"= {final_prob.numerator} / {final_prob.denominator}")

print("\nSo, the final probability is:")
print(f"{final_prob.numerator}/{final_prob.denominator}")
