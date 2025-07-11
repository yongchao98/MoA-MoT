import numpy as np

# Define probabilities and odds
p_true = np.array([1/2, 1/4, 1/4])
b = np.array([4, 3, 3])

# Calculate W_star
f_A_star = (p_true[0]*(b[0]+1)-1)/b[0] # 3/8
w_star_val = p_true[0]*np.log(1+f_A_star*b[0]) + (p_true[1]+p_true[2])*np.log(1-f_A_star)

# Calculate W
f_A_actual = 7/44
f_B_actual = 17/44
F_actual = f_A_actual + f_B_actual

R_A = 1 - F_actual + f_A_actual * (b[0] + 1)
R_B = 1 - F_actual + f_B_actual * (b[1] + 1)
R_C = 1 - F_actual

w_actual_val = p_true[0] * np.log(R_A) + p_true[1] * np.log(R_B) + p_true[2] * np.log(R_C)

# --- Output the step-by-step calculation of the difference ---
print("The optimal growth rate W* is calculated from the true probabilities and optimal fractions:")
print("W* = (1/2) * log(1 + (3/8)*4) + (1/4)*log(1-3/8) + (1/4)*log(1-3/8)")
print(f"W* = log(5/4) = {w_star_val:.4f}\n")


print("The actual growth rate W is calculated from the mistaken fractions but true probabilities:")
print("W = (1/2) * log(5/4) + (1/4) * log(2) + (1/4) * log(5/11)")
print(f"W = {w_actual_val:.4f}\n")

print("The difference W* - W is:")
# Final expression with all numbers shown
print("W* - W = [log(5/4)] - [(1/2)*log(5/4) + (1/4)*log(2) + (1/4)*log(5/11)]")
# Simplification steps
print("       = (1/2)*log(5/4) - (1/4)*log(2) - (1/4)*log(5/11)")
print("       = (1/4) * [2*log(5/4) - log(2) - log(5/11)]")
print("       = (1/4) * [log((5/4)^2) - log(2 * 5/11)]")
print("       = (1/4) * [log(25/16) - log(10/11)]")
print("       = (1/4) * log((25/16) / (10/11))")
print("       = (1/4) * log(275/160)")
print("       = (1/4) * log(55/32)\n")

final_value = w_star_val - w_actual_val
print(f"The final numerical value is:")
print(f"{final_value:.6f}")
<<<0.135394>>>