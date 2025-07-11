import math

# Define the probability of Theo winning a single game, derived as 1/6.
p_theo_wins_num = 1
p_theo_wins_den = 6

# The probability of Theo not winning is 1 - p_T = 5/6.
p_theo_not_wins_num = p_theo_wins_den - p_theo_wins_num
p_theo_not_wins_den = p_theo_wins_den

# The number of games Theo must not win in a row.
num_games = 4

print(f"The probability of Theo winning a single game is {p_theo_wins_num}/{p_theo_wins_den}.")
print(f"Therefore, the probability of Theo NOT winning a single game is {p_theo_not_wins_num}/{p_theo_not_wins_den}.")
print("")
print("The problem asks for the probability that Theo wins for the first time only after at least five games.")
print("This is the same as the probability that Theo does not win in the first four games.")
print(f"The calculation is ({p_theo_not_wins_num}/{p_theo_not_wins_den}) ^ {num_games}.")
print("")

# Calculate the numerator and denominator of the final probability.
final_numerator = p_theo_not_wins_num ** num_games
final_denominator = p_theo_not_wins_den ** num_games

print("The final probability is calculated as:")
print(f"{p_theo_not_wins_num}^{num_games} / {p_theo_not_wins_den}^{num_games} = {final_numerator} / {final_denominator}")

# Also print the decimal representation.
decimal_probability = final_numerator / final_denominator
print(f"\nAs a decimal, this is approximately: {decimal_probability}")
