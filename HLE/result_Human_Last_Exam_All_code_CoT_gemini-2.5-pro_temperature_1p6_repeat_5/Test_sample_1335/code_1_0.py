import sys

# Parameters from the problem statement
tax_on_capital_X_pct = 2.0
change_consumption_Y_pct = 3.0
labour_share_X = 2/3
own_price_elasticity_X = -2
own_price_elasticity_Y = -1
income_elasticity_X = 1
income_elasticity_Y = 1

# Convert percentages to decimals for calculation
tax_on_capital_X = tax_on_capital_X_pct / 100
C_hat_Y = change_consumption_Y_pct / 100

# --- Step 1: Percentage change in nominal wage (w_hat) ---
# In a small open economy with fixed P_Y and r, the wage 'w' is determined by the
# zero-profit condition in sector Y. Thus, w is fixed.
w_hat = 0.0

# --- Step 2: Percentage change in price of good X (P_hat_X) ---
# Using the price equation for sector X: P_hat_X = (labour_share_X * w_hat) + (capital_share_X * r_hat_X)
capital_share_X = 1 - labour_share_X
r_hat_X = tax_on_capital_X
P_hat_X = labour_share_X * w_hat + capital_share_X * r_hat_X

# --- Step 3: Percentage change in nominal income (I_hat) ---
# Start with the demand equation for Y: C_hat_Y = (eps_YX * P_hat_X) + (eps_YY * P_hat_Y) + (eta_Y * I_hat)
# From the homogeneity property of demand: eps_YX + eps_YY + eta_Y = 0
eps_YX = -own_price_elasticity_Y - income_elasticity_Y
P_hat_Y = 0.0  # Price of good Y is fixed (world price)
# Rearranging to solve for I_hat:
I_hat = (C_hat_Y - eps_YX * P_hat_X - own_price_elasticity_Y * P_hat_Y) / income_elasticity_Y

# --- Step 4: Percentage change in consumption of good X (C_hat_X) ---
# Use the demand equation for X: C_hat_X = (eps_XX * P_hat_X) + (eps_XY * P_hat_Y) + (eta_X * I_hat)
# From the homogeneity property of demand: eps_XX + eps_XY + eta_X = 0
eps_XY = -own_price_elasticity_X - income_elasticity_X
C_hat_X = own_price_elasticity_X * P_hat_X + eps_XY * P_hat_Y + income_elasticity_X * I_hat

# --- Final Output ---
# Convert results to percentages
w_hat_pct = w_hat * 100
P_hat_X_pct = P_hat_X * 100
C_hat_X_pct = C_hat_X * 100

# Print the key calculations as requested
print("Calculation for the percentage change in the price of X (P_hat_X):")
print(f"P_hat_X % = (labour_share_X * w_hat %) + (capital_share_X * r_hat_X %)")
print(f"P_hat_X % = ({labour_share_X:.4f} * {w_hat_pct:.4f} %) + ({capital_share_X:.4f} * {tax_on_capital_X_pct:.4f} %) = {P_hat_X_pct:.4f} %")
print("-" * 20)
print("Calculation for the percentage change in consumption of X (C_hat_X):")
print(f"C_hat_X % = (own_price_elasticity_X * P_hat_X %) + (income_elasticity_X * I_hat %)")
print(f"C_hat_X % = ({own_price_elasticity_X} * {P_hat_X_pct:.4f} %) + ({income_elasticity_X} * {I_hat * 100:.4f} %) = {C_hat_X_pct:.4f} %")
print("-" * 20)

# Final answer in comma-separated format
final_answer_string = f"{w_hat_pct},{P_hat_X_pct},{C_hat_X_pct}"
print("The percentage changes for nominal wage, price of good X, and consumption of good X are:")
print(final_answer_string)

# The final answer format for the platform
sys.stdout.write(f'<<<{final_answer_string}>>>')