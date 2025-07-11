import math

# Step 1: Define the game state
n1 = 136 # Number of 1-euro coins
n2 = 87  # Number of 2-euro coins

# Total number of coins
total_coins = n1 + n2

# Total value of all coins
total_value = n1 * 1 + n2 * 2

# Step 2: Determine player turns
# Since total_coins is odd, Player 1 gets one more turn than Player 2.
p1_coins_count = math.ceil(total_coins / 2)
p2_coins_count = math.floor(total_coins / 2)

print(f"There are {n1} 1-euro coins and {n2} 2-euro coins.")
print(f"The total number of coins is {total_coins}, which is an odd number.")
print(f"Player 1 will pick {p1_coins_count} coins in total.")
print(f"Player 2 will pick {p2_coins_count} coins in total.")
print("-" * 30)

# Step 3: Formulate the Winning Condition
# Let n2_1 be the number of 2-euro coins Player 1 collects.
# Let n1_1 be the number of 1-euro coins Player 1 collects.
# Let n2_2 be the number of 2-euro coins Player 2 collects.
# Let n1_2 be the number of 1-euro coins Player 2 collects.
#
# We know the following relationships:
# n1_1 + n2_1 = p1_coins_count  (Total coins for Player 1)
# n1_2 + n2_2 = p2_coins_count  (Total coins for Player 2)
# n1_1 + n1_2 = n1            (Total 1-euro coins)
# n2_1 + n2_2 = n2            (Total 2-euro coins)
#
# We can express Player 1's and Player 2's scores (S1, S2) in terms of n2_1.
# S1 = n1_1 * 1 + n2_1 * 2
#    = (p1_coins_count - n2_1) * 1 + n2_1 * 2
#    = p1_coins_count + n2_1
#
# S2 = Total Value - S1
#    = total_value - (p1_coins_count + n2_1)
#
# Player 1 wins if S1 > S2.
# p1_coins_count + n2_1 > total_value - (p1_coins_count + n2_1)
# 2 * (p1_coins_count + n2_1) > total_value
# 2 * p1_coins_count + 2 * n2_1 > total_value
# 2 * n2_1 > total_value - 2 * p1_coins_count

win_threshold_numerator = total_value - 2 * p1_coins_count
win_threshold = win_threshold_numerator / 2

print("Let n2_1 be the number of 2-euro coins collected by Player 1.")
print("The score for Player 1 (S1) can be expressed as: S1 = " + str(p1_coins_count) + " + n2_1")
print("The score for Player 2 (S2) can be expressed as: S2 = " + str(total_value) + " - S1 = " + str(total_value - p1_coins_count) + " - n2_1")
print("\nPlayer 1 wins if S1 > S2.")
print(f"This translates to the inequality: {p1_coins_count} + n2_1 > {total_value - p1_coins_count} - n2_1")
print(f"Simplifying the inequality, we get: 2 * n2_1 > {total_value - 2*p1_coins_count}")
print(f"Which finally means: n2_1 > {win_threshold}")
print("-" * 30)

# Step 4 & 5: Analysis and Conclusion
print(f"So, Player 1 wins if they collect more than {win_threshold} of the 2-euro coins.")
print(f"Since the number of coins must be an integer, Player 1 must collect at least {math.floor(win_threshold) + 1} of the {n2} available 2-euro coins to win.")
print(f"The game boils down to a fight for the majority of the {n2} 2-euro coins.")
print("Since there is an odd number of 2-euro coins (87), one player must collect more than the other. A tie is impossible.")
print("\nIn this game, the first player has 112 moves while the second has 111. Having an extra move, which is both the first and the last move of the game, provides a significant strategic advantage.")
print("This advantage allows the first player to better control the flow of the game and ensure they can collect the majority of the crucial 2-euro coins, regardless of the initial random arrangement.")
print("Therefore, it is always better to be the 1st player.")
