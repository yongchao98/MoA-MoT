import math

# --- Game Parameters ---
n1 = 136  # Number of 1-euro coins
n2 = 87   # Number of 2-euro coins

# --- Positional Parameters ---
total_coins = n1 + n2
odd_positions = (total_coins + 1) // 2
even_positions = total_coins - odd_positions
total_value = (n1 * 1) + (n2 * 2)
win_threshold = total_value / 2

# --- Calculate the most likely scenario ---
# N: population size (total_coins)
# K: successes in population (n2, the 2-euro coins)
# n: sample size (odd_positions)
mode_o2 = math.floor(((odd_positions + 1) * (n2 + 1)) / (total_coins + 2))

# Let o2 be the number of 2-euro coins in odd positions.
# The most likely value for o2 is its mode.
o2 = mode_o2

# --- Calculate S_odd and S_even for the most likely scenario ---
# S_odd = (num_1euro_odd * 1) + (num_2euro_odd * 2)
# num_1euro_odd = odd_positions - o2
# num_2euro_odd = o2
s_odd = (odd_positions - o2) * 1 + o2 * 2
s_even = total_value - s_odd

print(f"Total coins: {total_coins}")
print(f"Coins for Player 1: {odd_positions}")
print(f"Coins for Player 2: {even_positions}")
print(f"Total value: {total_value} euros")
print(f"Value needed to win: >{win_threshold} euros\n")

print("--- Analysis of the Most Likely Scenario ---")
print(f"The most likely number of 2-euro coins in the {odd_positions} odd positions is {o2}.")

print("\nFor this scenario, we calculate the value of coins in odd vs. even positions:")
print(f"S_odd = (Number of odd positions - o2) * 1 + o2 * 2")
print(f"S_odd = ({odd_positions} - {o2}) * 1 + {o2} * 2 = {s_odd}")
print(f"S_even = Total Value - S_odd = {total_value} - {s_odd} = {s_even}\n")

print("Player 2, being rational, sees S_odd > S_even and will play to claim the odd-position coins.")
print("Player 1's score is then dependent on their first move (value v_first).")
print("P1_score = S_even + v_first")
print("P2_score = S_odd - v_first\n")

# Case 1: Player 1 can pick a 2-euro coin from an end.
v_first_case1 = 2
p1_score_case1 = s_even + v_first_case1
p2_score_case1 = s_odd - v_first_case1
print(f"If P1's first pick has value {v_first_case1}:")
print(f"P1_score = {s_even} + {v_first_case1} = {p1_score_case1} (P1 Wins)")
print(f"P2_score = {s_odd} - {v_first_case1} = {p2_score_case1}")

# Case 2: Player 1 must pick a 1-euro coin from the ends.
v_first_case2 = 1
p1_score_case2 = s_even + v_first_case2
p2_score_case2 = s_odd - v_first_case2
print(f"\nIf P1's first pick has value {v_first_case2}:")
print(f"P1_score = {s_even} + {v_first_case2} = {p1_score_case2} (Draw)")
print(f"P2_score = {s_odd} - {v_first_case2} = {p2_score_case2}\n")

print("Conclusion: In the most likely arrangement, Player 1 can always secure at least a draw,")
print("and wins if a 2-euro coin is available at either end. This makes Player 1 the better role.")
