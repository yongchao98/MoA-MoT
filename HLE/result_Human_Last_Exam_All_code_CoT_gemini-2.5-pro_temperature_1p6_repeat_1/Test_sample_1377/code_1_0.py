import math

# --- Part 1: Optimal Growth Rate (W*) ---

# True probabilities
p_A, p_B, p_C = 1/2, 1/4, 1/4

# Payout odds
b_A, b_B, b_C = 4, 3, 3

# Optimal fraction for A (only positive edge)
f_star_A = p_A - (1 - p_A) / b_A

# Return if A wins
R_star_A_wins = 1 + f_star_A * b_A
# Return if A loses
R_star_A_loses = 1 - f_star_A

# Optimal Growth Rate W*
W_star = p_A * math.log(R_star_A_wins) + (p_B + p_C) * math.log(R_star_A_loses)

print("--- Optimal Strategy (W*) ---")
print(f"Optimal fraction to bet on A: f_A* = {f_star_A:.4f}")
print(f"Equation for W*: W* = {p_A}*ln({R_star_A_wins:.4f}) + {p_B+p_C}*ln({R_star_A_loses:.4f})")
# The equation simplifies to ln(5/4)
W_star_simplified_term = 5/4
print(f"Simplified W* = ln({W_star_simplified_term})")
print(f"W* = {W_star:.5f}\n")


# --- Part 2: Actual Growth Rate (W) ---

# Fractions bet under mistaken beliefs
f_A = 7/44
f_B = 17/44

# Returns based on mistaken bets
R_A = 1 - f_B + f_A * b_A
R_B = 1 - f_A + f_B * b_B
R_C = 1 - f_A - f_B

# Actual growth rate W, using true probabilities
W = p_A * math.log(R_A) + p_B * math.log(R_B) + p_C * math.log(R_C)

print("--- Actual Strategy (W) ---")
print(f"Mistaken fraction on A: f_A = {f_A:.4f}")
print(f"Mistaken fraction on B: f_B = {f_B:.4f}")
print(f"Equation for W: W = ({p_A})*ln({R_A:.4f}) + ({p_B})*ln({R_B:.4f}) + ({p_C})*ln({R_C:.4f})")
print(f"W = {W:.5f}\n")


# --- Part 3: Difference (W* - W) ---

difference = W_star - W

# Final equation components
final_term1 = W_star_simplified_term
final_term2_1_val = R_A
final_term2_2_val = R_B
final_term2_3_val = R_C
final_term_pA = p_A
final_term_pB = p_B
final_term_pC = p_C
final_simplified_term = 55/32

print("--- Final Result ---")
print("Equation for W* - W:")
print(f"  ln({final_term1}) - [({final_term_pA})*ln({final_term2_1_val:.4f}) + ({final_term_pB})*ln({final_term2_2_val:.4f}) + ({final_term_pC})*ln({final_term2_3_val:.4f})]")
print(f"= {W_star:.5f} - {W:.5f}")
print(f"= {difference:.5f}")
print(f"\nSimplified symbolic answer: (1/4)*ln({final_simplified_term:.4f})")
print("The calculated difference W* - W is:")
print(difference)