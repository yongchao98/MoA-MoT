import math

# Step 1: Define the probability of Theo winning a single game.
# As derived from the random walk model, p = 1/6.
p_theo_wins = 1/6

# Step 2: Define the probability of Theo not winning a single game.
p_theo_not_wins = 1 - p_theo_wins

# Step 3: The event "Theo wins for the first time only after at least five games"
# is equivalent to "Theo does not win in the first four games".
num_games_without_win = 4

# Step 4: Calculate the probability of this event.
# Since the games are independent, we raise the single-game probability
# of not winning to the power of the number of games.
numerator = int(5**num_games_without_win)
denominator = int(6**num_games_without_win)
final_probability = p_theo_not_wins ** num_games_without_win

# Step 5: Print the result, showing the components of the equation.
print("The probability of Theo winning a single game is 1/6.")
print("The probability of Theo not winning a single game is 5/6.")
print("The probability of Theo not winning the first four games is (5/6)^4.")
print(f"Calculation: ({5}/{6})^{num_games_without_win} = {numerator}/{denominator}")
print(f"The final probability is {final_probability:.6f}")
