import math

# Step 1: Define the probabilities for a single game.
# Through solving a random walk model, it can be shown that:
# P(Alex wins) = 1/6
# P(Theo wins) = 1/6
# P(Draw) = 2/3
# We are interested in the probability of Theo winning.
p_theo_wins = 1/6

# Step 2: Define the event of interest.
# The event "Theo wins for the first time only after at least five games"
# is the same as "Theo does not win in game 1 AND does not win in game 2
# AND does not win in game 3 AND does not win in game 4".

# Step 3: Calculate the probability of the event.
# The probability that Theo does NOT win a single game is 1 - p_theo_wins.
p_not_theo_wins_numerator = 5
p_not_theo_wins_denominator = 6

# The number of games in the sequence where Theo must not win is 4.
num_games_must_not_win = 4

# The total probability is (5/6)^4.
final_numerator = p_not_theo_wins_numerator ** num_games_must_not_win
final_denominator = p_not_theo_wins_denominator ** num_games_must_not_win

# Calculate the final result.
final_probability = final_numerator / final_denominator

# Step 4: Print the results, showing the components of the final equation.
print(f"The probability of Theo winning a single game is 1/6.")
print(f"Therefore, the probability of Theo NOT winning a single game is {p_not_theo_wins_numerator}/{p_not_theo_wins_denominator}.")
print(f"We need to find the probability that Theo does not win in the first {num_games_must_not_win} games.")
print(f"The final equation is ({p_not_theo_wins_numerator}/{p_not_theo_wins_denominator}) ** {num_games_must_not_win}")
print(f"Which is calculated as: {final_numerator} / {final_denominator}")
print(f"The final probability is: {final_probability}")
