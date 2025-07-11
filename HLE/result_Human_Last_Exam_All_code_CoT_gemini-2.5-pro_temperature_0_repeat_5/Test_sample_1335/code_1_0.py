import sys
import io

# Redirect stdout to capture print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Given Parameters ---
# Cost Shares
theta_LX = 2/3  # Labor's share in X
theta_KX = 1/3  # Capital's share in X
theta_LY = 1/3  # Labor's share in Y
theta_KY = 2/3  # Capital's share in Y

# The Shock
hat_r_X = 0.02  # 2% tax on capital return in sector X

# Demand Elasticities
eta_XX = -2  # Own-price elasticity for X
eta_YY = -1  # Own-price elasticity for Y
eta_I_X = 1   # Income elasticity for X
eta_I_Y = 1   # Income elasticity for Y

# Given Outcome
hat_C_Y = 0.03 # 3% increase in consumption of Y

# --- Model Assumptions ---
# Small Open Economy implies P_Y and r are fixed
hat_P_Y = 0
hat_r_Y = 0 # No tax in sector Y, world r is fixed

# --- Step 1: Calculate Percentage Change in Nominal Wage (hat_w) ---
print("Step 1: Calculate the percentage change in nominal wage (hat_w)")
print("The zero-profit condition in sector Y is: hat(P_Y) = theta_LY * hat(w) + theta_KY * hat(r_Y)")
print(f"Plugging in the values: {hat_P_Y} = {theta_LY:.2f} * hat(w) + {theta_KY:.2f} * {hat_r_Y}")
# 0 = (1/3) * hat_w + (2/3) * 0
hat_w = 0
print(f"This simplifies to hat(w) = {hat_w:.4f}")
print("-" * 20)

# --- Step 2: Calculate Percentage Change in Price of Good X (hat_P_X) ---
print("Step 2: Calculate the percentage change in the price of good X (hat_P_X)")
print("The zero-profit condition in sector X is: hat(P_X) = theta_LX * hat(w) + theta_KX * hat(r_X)")
# hat_P_X = (2/3) * 0 + (1/3) * 0.02
hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
print(f"Plugging in the values: hat(P_X) = {theta_LX:.2f} * {hat_w:.4f} + {theta_KX:.2f} * {hat_r_X}")
print(f"This gives hat(P_X) = {hat_P_X:.4f}")
print("-" * 20)

# --- Step 3: Calculate Percentage Change in Consumption of Good X (hat_C_X) ---
print("Step 3: Calculate the percentage change in the consumption of good X (hat_C_X)")

# First, find the percentage change in income (hat_I) using the demand for Y
print("First, we find the change in income (hat_I) using the information for good Y.")
print("The demand function for Y is: hat(C_Y) = eta_YX * hat(P_X) + eta_YY * hat(P_Y) + eta_I_Y * hat(I)")
# From homogeneity of demand: eta_YX + eta_YY + eta_I_Y = 0
eta_YX = -eta_YY - eta_I_Y
print(f"From homogeneity, eta_YX = -eta_YY - eta_I_Y = -({eta_YY}) - {eta_I_Y} = {eta_YX}")
# 0.03 = 0 * hat_P_X + (-1) * 0 + 1 * hat_I
hat_I = (hat_C_Y - eta_YX * hat_P_X - eta_YY * hat_P_Y) / eta_I_Y
print(f"Plugging into the demand function: {hat_C_Y} = {eta_YX} * {hat_P_X:.4f} + {eta_YY} * {hat_P_Y} + {eta_I_Y} * hat(I)")
print(f"Solving for hat(I) gives: hat(I) = {hat_I:.4f}")
print("")

# Now, find the percentage change in consumption of X (hat_C_X)
print("Now, we can find the change in consumption of X (hat_C_X).")
print("The demand function for X is: hat(C_X) = eta_XX * hat(P_X) + eta_XY * hat(P_Y) + eta_I_X * hat(I)")
# From homogeneity of demand: eta_XX + eta_XY + eta_I_X = 0
eta_XY = -eta_XX - eta_I_X
print(f"From homogeneity, eta_XY = -eta_XX - eta_I_X = -({eta_XX}) - {eta_I_X} = {eta_XY}")
# hat_C_X = (-2) * hat_P_X + 1 * 0 + 1 * hat_I
hat_C_X = eta_XX * hat_P_X + eta_XY * hat_P_Y + eta_I_X * hat_I
print(f"Plugging into the demand function: hat(C_X) = {eta_XX} * {hat_P_X:.4f} + {eta_XY} * {hat_P_Y} + {eta_I_X} * {hat_I:.4f}")
print(f"Solving for hat(C_X) gives: hat(C_X) = {hat_C_X:.4f}")
print("-" * 20)

# --- Final Answer ---
# Convert proportional changes to percentage changes
hat_w_percent = hat_w * 100
hat_P_X_percent = hat_P_X * 100
hat_C_X_percent = hat_C_X * 100

# Restore stdout
sys.stdout = old_stdout
# Print the captured output to the console
print(captured_output.getvalue())

# Print the final comma-separated result
print("Final comma-separated percentage changes (wage, price of X, consumption of X):")
final_answer_string = f"{hat_w_percent},{hat_P_X_percent},{hat_C_X_percent}"
print(final_answer_string)
print(f"<<<{final_answer_string}>>>")