import sys
from io import StringIO

# Redirect stdout to capture the print output
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

# --- Step 1: Define Parameters from the problem ---
# Cost share of labor in Y
theta_LY = 1/3
# Cost share of capital in X
theta_KX = 1/3
# Tax on capital in X, which is the percentage change in the rental rate for sector X
hat_r_X = 0.02
# Given percentage change in consumption of Y
hat_C_Y = 0.03
# Own-price elasticity of demand for Y
epsilon_YY = -1
# Income elasticity of demand for Y
eta_Y = 1
# Own-price elasticity of demand for X
epsilon_XX = -2
# Income elasticity of demand for X
eta_X = 1
# Percentage change in price of Y is 0 (freely traded)
hat_P_Y = 0
# Rental rate change in Y is 0 (tax only on X)
hat_r_Y = 0

# --- Step 2: Calculate the percentage change in nominal wage (hat_w) ---
# The zero-profit condition for sector Y is: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
# Since hat_P_Y = 0 and hat_r_Y = 0, the equation becomes 0 = (1/3) * hat_w + (2/3) * 0
hat_w = 0
print("Equation for percentage change in nominal wage (w):")
print(f"{hat_w} = 0 (derived from 0 = {theta_LY:.2f}*w + (1-{theta_LY:.2f})*0)")
print("-" * 20)

# --- Step 3: Calculate the percentage change in the price of good X (hat_P_X) ---
# The zero-profit condition for sector X is: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
# Since hat_w = 0, the equation becomes hat_P_X = theta_KX * hat_r_X
hat_P_X = theta_KX * hat_r_X
print("Equation for percentage change in price of good X (P_X):")
print(f"{hat_P_X:.4f} = (1 - {2/3:.2f}) * {hat_w} + {theta_KX:.2f} * {hat_r_X}")
print("-" * 20)

# --- Step 4: Calculate the percentage change in nominal income (hat_I) ---
# From homogeneity of demand, the sum of price and income elasticities is zero.
# For good Y: epsilon_YX + epsilon_YY + eta_Y = 0
# This gives epsilon_YX = -epsilon_YY - eta_Y
epsilon_YX = -epsilon_YY - eta_Y
# The demand function for Y is: hat_C_Y = epsilon_YX * hat_P_X + epsilon_YY * hat_P_Y + eta_Y * hat_I
# We can solve for hat_I: hat_I = (hat_C_Y - epsilon_YX * hat_P_X - epsilon_YY * hat_P_Y) / eta_Y
hat_I = (hat_C_Y - epsilon_YX * hat_P_X - epsilon_YY * hat_P_Y) / eta_Y
print("Equation for percentage change in nominal income (I):")
print(f"{hat_I:.4f} = ({hat_C_Y} - {epsilon_YX} * {hat_P_X:.4f}) / {eta_Y}")
print("-" * 20)

# --- Step 5: Calculate the percentage change in consumption of good X (hat_C_X) ---
# The demand function for X is: hat_C_X = epsilon_XX * hat_P_X + epsilon_XY * hat_P_Y + eta_X * hat_I
# Since hat_P_Y = 0, the equation is: hat_C_X = epsilon_XX * hat_P_X + eta_X * hat_I
hat_C_X = epsilon_XX * hat_P_X + eta_X * hat_I
print("Equation for percentage change in consumption of good X (C_X):")
print(f"{hat_C_X:.4f} = {epsilon_XX} * {hat_P_X:.4f} + {eta_X} * {hat_I:.4f}")
print("-" * 20)

# --- Final Answer ---
# The final comma-separated values for (hat_w, hat_P_X, hat_C_X)
final_answer_string = f"{hat_w},{hat_P_X},{hat_C_X}"
print("Final comma-separated values (wage, price of X, consumption of X):")
print(final_answer_string)

# Restore stdout
sys.stdout = old_stdout
# Print the captured output to the actual console
output = captured_output.getvalue()
print(output)
final_answer = output.strip().split('\n')[-1]
print(f"<<<{final_answer}>>>")