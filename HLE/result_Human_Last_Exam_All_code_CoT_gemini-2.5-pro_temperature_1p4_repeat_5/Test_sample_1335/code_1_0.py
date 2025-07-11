import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# --- Model Parameters ---
# Cost shares
theta_LX = 2/3  # Labor's share of cost in sector X
theta_KX = 1 - theta_LX  # Capital's share of cost in sector X
theta_LY = 1/3  # Labor's share of cost in sector Y
theta_KY = 1 - theta_LY  # Capital's share of cost in sector Y

# Demand Elasticities
eta_pX = -2.0  # Own price elasticity of demand for good X
eta_pY = -1.0  # Own price elasticity of demand for good Y
eta_IX = 1.0   # Income elasticity of demand for good X
eta_IY = 1.0   # Income elasticity of demand for good Y

# Shocks and Known Outcomes
tax_rate = 0.02         # 2% tax on capital in sector X
hat_r_X = tax_rate      # Percentage change in capital cost in X
hat_C_Y = 0.03          # 3% increase in consumption of Y

# Fixed world prices
hat_P_Y = 0.0  # Percentage change in price of good Y is 0
hat_r_Y = 0.0  # Percentage change in world return to capital is 0

# --- Step 1: Calculate the percentage change in nominal wage (hat_w) ---
# From the zero-profit condition in sector Y: hat(P_Y) = theta_LY * hat(w) + theta_KY * hat(r_Y)
# 0 = (1/3) * hat(w) + (2/3) * 0
# This implies hat(w) = 0
hat_w = 0.0
print("Step 1: Calculate the percentage change in nominal wage (hat_w)")
print(f"Equation: hat(P_Y) = theta_LY * hat(w) + theta_KY * hat(r_Y)")
print(f"Plugging in values: {hat_P_Y} = {theta_LY:.4f} * hat(w) + {theta_KY:.4f} * {hat_r_Y}")
print(f"Result: hat(w) = {hat_w:.4f}\n")


# --- Step 2: Calculate the percentage change in the price of good X (hat_P_X) ---
# From the zero-profit condition in sector X: hat(P_X) = theta_LX * hat(w) + theta_KX * hat(r_X)
hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
print("Step 2: Calculate the percentage change in the price of good X (hat_P_X)")
print(f"Equation: hat(P_X) = theta_LX * hat(w) + theta_KX * hat(r_X)")
print(f"Plugging in values: hat(P_X) = {theta_LX:.4f} * {hat_w:.4f} + {theta_KX:.4f} * {hat_r_X}")
print(f"Result: hat(P_X) = {hat_P_X:.4f}\n")


# --- Step 3: Calculate the percentage change in income (hat_I) ---
# From the demand function for good Y: hat(C_Y) = eta_pY * hat(P_Y) + eta_IY * hat(I)
# hat(I) = (hat(C_Y) - eta_pY * hat(P_Y)) / eta_IY
hat_I = (hat_C_Y - eta_pY * hat_P_Y) / eta_IY
print("Step 3: Calculate the percentage change in income (hat_I)")
print(f"Equation: hat(C_Y) = eta_pY * hat(P_Y) + eta_IY * hat(I)")
print(f"Plugging in values: {hat_C_Y} = {eta_pY} * {hat_P_Y} + {eta_IY} * hat(I)")
print(f"Result: hat(I) = {hat_I:.4f}\n")


# --- Step 4: Calculate the percentage change in consumption of good X (hat_C_X) ---
# From the demand function for good X: hat(C_X) = eta_pX * hat(P_X) + eta_IX * hat(I)
hat_C_X = eta_pX * hat_P_X + eta_IX * hat_I
print("Step 4: Calculate the percentage change in consumption of good X (hat_C_X)")
print(f"Equation: hat(C_X) = eta_pX * hat(P_X) + eta_IX * hat(I)")
print(f"Plugging in values: hat(C_X) = {eta_pX} * {hat_P_X:.4f} + {eta_IX} * {hat_I:.4f}")
print(f"Result: hat(C_X) = {hat_C_X:.4f}\n")

# --- Final Answer ---
# Convert proportional changes to percentages for the final output
wage_change_pct = hat_w * 100
price_X_change_pct = hat_P_X * 100
consumption_X_change_pct = hat_C_X * 100

print("Final Answer: The percentage changes are as follows:")
print(f"Percentage change in nominal wage: {wage_change_pct:.4f}%")
print(f"Percentage change in price of good X: {price_X_change_pct:.4f}%")
print(f"Percentage change in consumption of good X: {consumption_X_change_pct:.4f}%")

# Restore the original stdout
sys.stdout = original_stdout
# Get the content of the buffer
output = buffer.getvalue()
# Print the buffer content
print(output)

# Print the final comma-separated values as requested by the user
final_answer = f"{wage_change_pct},{price_X_change_pct},{consumption_X_change_pct}"
print(final_answer)
print(f"<<<{final_answer}>>>")