import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Given Parameters ---
# Elasticities
eta_I = 1.0  # Income elasticity for both goods
eta_XX = -2.0  # Own-price elasticity for X
eta_YY = -1.0  # Own-price elasticity for Y

# Cost Shares
theta_LX = 2/3  # Labor's share of cost in sector X
theta_KX = 1 - theta_LX  # Capital's share of cost in sector X
theta_LY = 1/3  # Labor's share of cost in sector Y
theta_KY = 1 - theta_LY  # Capital's share of cost in sector Y

# Shocks and results
hat_rX = 0.02  # 2% tax on capital in sector X
hat_CY = 0.03  # 3% increase in consumption of Y

# --- Step 1: Calculate percentage change in wage (hat_w) ---
# From the zero-profit condition in the traded sector Y:
# hat_PY = theta_LY * hat_w + theta_KY * hat_rY
# Since Y is traded and K is mobile, P_Y and r are fixed.
# So, hat_PY = 0 and hat_rY = 0.
# 0 = theta_LY * hat_w + theta_KY * 0
# This implies hat_w = 0.
hat_w = 0.0

# --- Step 2: Calculate percentage change in price of X (hat_PX) ---
# From the zero-profit condition in sector X:
# hat_PX = theta_LX * hat_w + theta_KX * hat_rX
hat_PX = theta_LX * hat_w + theta_KX * hat_rX

# --- Step 3: Calculate percentage change in consumption of X (hat_CX) ---
# This requires finding the change in income, hat_I.
# We use the demand equations and demand theory.

# First, find the cross-price elasticity eta_YX using Cournot aggregation:
# s_X * eta_XX + s_Y * eta_YX = -s_X
# To solve this, we need the expenditure shares s_X and s_Y.
# As they are not given, we make the simplifying assumption of equal shares.
s_X = 0.5
s_Y = 0.5

# With s_X and s_Y, we can find eta_YX
# s_Y * eta_YX = -s_X - s_X * eta_XX
# eta_YX = (-s_X - s_X * eta_XX) / s_Y
eta_YX = (-s_X * (1 + eta_XX)) / s_Y

# Now use the demand equation for Y to find hat_I:
# hat_CY = eta_YX * hat_PX + eta_YY * hat_PY + eta_I * hat_I
# We know hat_PY = 0.
# hat_CY = eta_YX * hat_PX + hat_I
# So, hat_I = hat_CY - eta_YX * hat_PX
hat_I = hat_CY - eta_YX * hat_PX

# Finally, use the demand equation for X to find hat_CX:
# hat_CX = eta_XX * hat_PX + eta_XY * hat_PY + eta_I * hat_I
# We know hat_PY = 0.
# hat_CX = eta_XX * hat_PX + hat_I
hat_CX = eta_XX * hat_PX + hat_I

# --- Output the results ---
# Convert changes to percentages for the final answer
hat_w_percent = hat_w * 100
hat_PX_percent = hat_PX * 100
hat_CX_percent = hat_CX * 100

# We need to print the calculation steps as requested by the prompt.
print(f"1. Percentage change in nominal wage (w):")
print(f"The zero profit condition for sector Y is: dPY/PY = theta_LY * (dw/w) + theta_KY * (dr_Y/r_Y)")
print(f"Since Y is traded and K is mobile, dPY/PY = 0 and dr_Y/r_Y = 0.")
print(f"0 = {theta_LY:.2f} * (dw/w) + {theta_KY:.2f} * 0")
print(f"Therefore, the percentage change in nominal wage is {hat_w_percent:.3f}%.")
print("-" * 20)

print(f"2. Percentage change in price of good X (PX):")
print(f"The zero profit condition for sector X is: dPX/PX = theta_LX * (dw/w) + theta_KX * (dr_X/r_X)")
print(f"We know dw/w = {hat_w:.3f} and the tax implies dr_X/r_X = {hat_rX:.2f}.")
print(f"dPX/PX = {theta_LX:.2f} * {hat_w:.3f} + {theta_KX:.2f} * {hat_rX:.2f}")
print(f"dPX/PX = {hat_PX:.6f}")
print(f"Therefore, the percentage change in the price of good X is {hat_PX_percent:.3f}%.")
print("-" * 20)

print(f"3. Percentage change in consumption of good X (CX):")
print(f"The demand function for Y is: dCY/CY = eta_YX * (dPX/PX) + eta_I * (dI/I)")
print(f"We are given dCY/CY = {hat_CY:.2f}. We need to find dI/I.")
print(f"Assuming equal budget shares (s_X = s_Y = 0.5), Cournot aggregation (s_X*eta_XX + s_Y*eta_YX = -s_X) gives eta_YX = {eta_YX:.2f}.")
print(f"{hat_CY:.2f} = {eta_YX:.2f} * {hat_PX:.6f} + 1 * (dI/I)")
hat_I_val = hat_CY - eta_YX * hat_PX
print(f"dI/I = {hat_CY:.2f} - {eta_YX * hat_PX:.6f} = {hat_I_val:.6f}")
print(f"The demand function for X is: dCX/CX = eta_XX * (dPX/PX) + eta_I * (dI/I)")
print(f"dCX/CX = {eta_XX:.2f} * {hat_PX:.6f} + 1 * {hat_I_val:.6f}")
print(f"dCX/CX = {eta_XX * hat_PX:.6f} + {hat_I_val:.6f} = {hat_CX:.6f}")
print(f"Therefore, the percentage change in the consumption of good X is {hat_CX_percent:.3f}%.")
print("-" * 20)

final_answer_str = f"{hat_w_percent:.3f},{hat_PX_percent:.3f},{hat_CX_percent:.3f}"
print(f"Final comma-separated answer: {final_answer_str}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the output to the actual console
print(output)
print(f'<<<{final_answer_str}>>>')