import math

# Step 1: Determine the probability of Theo winning a single game, P(T).
# The game can be modeled as a random walk. For Theo to win, the walk must start
# by going to state d=-1 (first toss is Tails, prob 1/2), and then from d=-1,
# it must reach the winning state d=-3 before reaching the draw state d=0.
# By symmetry with Alex's win condition, we can calculate P(A).
# For Alex to win from d=1 (after a Heads), he must reach d=3 before d=0.
# This is a Gambler's Ruin problem: P(win from 1) = (1-0)/(3-0) = 1/3.
# So, P(A) = P(first toss H) * P(win from 1) = (1/2) * (1/3) = 1/6.
p_A_num = 1
p_A_den = 6

# By symmetry, P(T) = P(A).
p_T_num = p_A_num
p_T_den = p_A_den

print(f"The probability of Theo winning a single game is P(T) = {p_T_num}/{p_T_den}.")

# Step 2: Find the probability of Theo NOT winning a game.
# P(Not T) = 1 - P(T)
p_not_T_num = p_T_den - p_T_num
p_not_T_den = p_T_den

print(f"The probability of Theo not winning a single game is {p_not_T_num}/{p_not_T_den}.")

# Step 3: Calculate the probability of this event occurring for 4 consecutive games.
# "Theo wins for the first time only after at least five games" means he does
# not win game 1, AND not win game 2, AND not win game 3, AND not win game 4.
n = 4

print(f"The problem asks for the probability that Theo doesn't win in the first {n} games.")

# The final probability is (P(Not T))^n
result_num = p_not_T_num ** n
result_den = p_not_T_den ** n
final_prob = result_num / result_den

print(f"The calculation is ({p_not_T_num}/{p_not_T_den})^{n} = {result_num}/{result_den}")
print(f"The probability that Theo wins for the first time only after at least five games is {final_prob:.7f}")
