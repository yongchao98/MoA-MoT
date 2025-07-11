import math

# Step 1: Define the known probabilities for a single game.
# Based on the random walk analysis, the probability of Theo winning is 1/6.
p_theo_wins_num = 1
p_theo_wins_den = 6

# Step 2: The problem asks for the probability that Theo's first win occurs
# on or after the 5th game. This is equivalent to the probability that
# Theo does not win the first four games.

# The probability of Theo NOT winning a single game is (1 - p).
p_theo_not_wins_num = p_theo_wins_den - p_theo_wins_num
p_theo_not_wins_den = p_theo_wins_den

# The number of games Theo must not win is 4.
n = 4

# Step 3: Calculate the final probability, which is (1 - p)^n.
final_prob_num = p_theo_not_wins_num ** n
final_prob_den = p_theo_not_wins_den ** n

# Step 4: Print the explanation and the final result.
print("The probability of Theo winning a single game is 1/6.")
print(f"The probability of Theo not winning a single game is 1 - {p_theo_wins_num}/{p_theo_wins_den} = {p_theo_not_wins_num}/{p_theo_not_wins_den}.")
print("\nThe probability that Theo's first win occurs after at least five games is the probability that he does not win the first four games.")
print(f"This is calculated as ({p_theo_not_wins_num}/{p_theo_not_wins_den})^{n}.")

print("\nFinal Equation:")
print(f"({p_theo_not_wins_num}/{p_theo_not_wins_den})^{n} = ({p_theo_not_wins_num}^{n}) / ({p_theo_not_wins_den}^{n}) = {final_prob_num}/{final_prob_den}")

# Calculate and print the decimal value
final_prob_decimal = final_prob_num / final_prob_den
print(f"\nAs a decimal, the probability is approximately: {final_prob_decimal:.5f}")
