import math

# Initial setup
num_ones = 136
num_twos = 87

# This specific arrangement has all 2-euro coins first, then all 1-euro coins.
coins = [2] * num_twos + [1] * num_ones

# Pointers to the left and right ends of the line
left = 0
right = len(coins) - 1

# Variables to store the number of each coin type collected by each player
p1_ones, p1_twos = 0, 0
p2_ones, p2_twos = 0, 0

# The game simulation
for turn in range(len(coins)):
    is_p1_turn = (turn % 2 == 0)

    # Optimal strategy for any player is to take the highest value coin available.
    # In this setup, the left coin is always >= the right coin for the first 87 moves.
    # After that, all coins are 1s, so the choice doesn't matter. We'll default to left.
    if coins[left] > coins[right]:
        picked_coin = coins[left]
        left += 1
    else:
        picked_coin = coins[left]
        left += 1

    # Add the collected coin to the respective player's count
    if is_p1_turn:
        if picked_coin == 1:
            p1_ones += 1
        else:
            p1_twos += 1
    else: # P2's turn
        if picked_coin == 1:
            p2_ones += 1
        else:
            p2_twos += 1

# Calculate the final scores
p1_score = p1_twos * 2 + p1_ones * 1
p2_score = p2_twos * 2 + p2_ones * 1

print("Demonstration with an arrangement where all 2-euro coins are on one side:")
print(f"Player 1 collects {p1_twos} 2-euro coins and {p1_ones} 1-euro coins.")
print(f"Player 2 collects {p2_twos} 2-euro coins and {p2_ones} 1-euro coins.")
print("")
print("Final Score Calculation:")
print(f"Player 1 Value = {p1_twos} * 2 + {p1_ones} * 1 = {p1_score}")
print(f"Player 2 Value = {p2_twos} * 2 + {p2_ones} * 1 = {p2_score}")
print("")
print("Result: Player 1 wins.")