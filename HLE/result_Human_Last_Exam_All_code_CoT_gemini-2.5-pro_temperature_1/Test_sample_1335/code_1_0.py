import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Given Parameters ---
# Elasticities of Substitution
sigma_X = -2  # Note: Non-standard value, but not needed for the solution.
sigma_Y = -1  # Note: Non-standard value, but not needed for the solution.

# Demand Elasticities
eta_IX = 1.0  # Income elasticity for X
eta_IY = 1.0  # Income elasticity for Y
eps_XX = -2.0 # Own-price elasticity for X
eps_YY = -1.0 # Own-price elasticity for Y

# Factor share allocations
lambda_LX = 3/4
lambda_KX = 1/4

# Cost shares
theta_LX = 2/3
theta_KX = 1 - theta_LX
theta_LY = 1/3
theta_KY = 1 - theta_LY

# Shock
tax_rate = 0.02
r_X_hat = tax_rate # Percentage change in return to capital in sector X

# Known outcome (not needed for the direct solution but given)
C_Y_hat_given = 0.03

# --- Step-by-step Calculation ---

# Step 1: Calculate the percentage change in nominal wage (w_hat)
# From the zero-profit condition in sector Y: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat
# In a small open economy, P_Y is fixed (P_Y_hat = 0) and the return to capital
# in the untaxed sector is fixed (r_Y_hat = 0).
# So, 0 = theta_LY * w_hat + theta_KY * 0
# This implies w_hat must be 0.
w_hat = 0.0
print(f"Calculation for percentage change in nominal wage (w_hat):")
print(f"From 0 = {theta_LY:.4f} * w_hat + {theta_KY:.4f} * 0, we get w_hat = {w_hat}")
print("-" * 20)

# Step 2: Calculate the percentage change in the price of good X (P_X_hat)
# From the zero-profit condition in sector X: P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
print(f"Calculation for percentage change in price of good X (P_X_hat):")
print(f"P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat")
print(f"P_X_hat = {theta_LX:.4f} * {w_hat} + {theta_KX:.4f} * {r_X_hat} = {P_X_hat:.8f}")
print("-" * 20)

# Step 3: Calculate the percentage change in consumption of good X (C_X_hat)
# From the demand function: C_X_hat = eps_XX * P_X_hat + eta_IX * I_hat
# We assume national income I = w*L (income from the fixed factor).
# So, I_hat = w_hat + L_hat. With fixed labor supply L_hat=0, so I_hat = w_hat.
I_hat = w_hat
C_X_hat = eps_XX * P_X_hat + eta_IX * I_hat
print(f"Calculation for percentage change in consumption of good X (C_X_hat):")
print(f"First, find income change (I_hat). Assuming I = w*L, I_hat = w_hat = {I_hat}")
print(f"Then, C_X_hat = eps_XX * P_X_hat + eta_IX * I_hat")
print(f"C_X_hat = {eps_XX} * {P_X_hat:.8f} + {eta_IX} * {I_hat} = {C_X_hat:.8f}")
print("-" * 20)

# --- Final Answer ---
# The results are the proportional changes. To get percentage changes, multiply by 100.
w_hat_percent = w_hat * 100
P_X_hat_percent = P_X_hat * 100
C_X_hat_percent = C_X_hat * 100

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
detailed_output = captured_output.getvalue()
print(detailed_output)

final_answer_string = f"{w_hat_percent},{P_X_hat_percent},{C_X_hat_percent}"
print(final_answer_string)
print(f"<<<{final_answer_string}>>>")