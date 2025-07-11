import math

# --- Given Parameters ---
tax_on_capital = 0.02 # 2% tax
pct_change_cy = 0.03   # Consumption of Y increases by 3%

# Cost shares
theta_lx = 2/3
theta_kx = 1/3
theta_ly = 1/3
theta_ky = 2/3

# Demand Elasticities
epsilon_xx = -2
epsilon_yy = -1
eta_ix = 1
eta_iy = 1

# --- Step-by-Step Calculation ---

# Step 1: Calculate the percentage change in nominal wage (w_hat)
# The zero-profit condition for the traded good Y is: P_y_hat = theta_ly * w_hat + theta_ky * r_hat
# As a small open economy, P_y_hat = 0 (fixed world price) and r_hat = 0 (fixed world capital return).
# So, 0 = (1/3) * w_hat + (2/3) * 0
w_hat = 0.0

# Step 2: Calculate the percentage change in the price of good X (px_hat)
# The zero-profit condition for good X is: Px_hat = theta_lx * w_hat + theta_kx * r_x_hat_firm
# The tax increases the firm's cost of capital in sector X by 2%, so r_x_hat_firm = 0.02.
# We know w_hat = 0.
# So, Px_hat = (2/3) * 0 + (1/3) * 0.02
px_hat = theta_kx * tax_on_capital

# Step 3: Calculate the percentage change in nominal income (i_hat)
# The demand function for Y is: Cy_hat = epsilon_yx * Px_hat + epsilon_yy * Py_hat + eta_iy * i_hat
# We know Py_hat = 0.
# From homogeneity of demand (epsilon_yx + epsilon_yy + eta_iy = 0), we find epsilon_yx.
# epsilon_yx + (-1) + 1 = 0  => epsilon_yx = 0
# The demand equation simplifies to: Cy_hat = eta_iy * i_hat
# 0.03 = 1 * i_hat
i_hat = pct_change_cy / eta_iy

# Step 4: Calculate the percentage change in consumption of good X (cx_hat)
# The demand function for X is: Cx_hat = epsilon_xx * Px_hat + epsilon_xy * Py_hat + eta_ix * i_hat
# We know Py_hat = 0.
# So, Cx_hat = epsilon_xx * Px_hat + eta_ix * i_hat
cx_hat = epsilon_xx * px_hat + eta_ix * i_hat

# --- Final Output ---
# Convert results to percentages for the final output
w_hat_pct = w_hat * 100
px_hat_pct = px_hat * 100
cx_hat_pct = cx_hat * 100

print("This script calculates the percentage changes in nominal wage, price of good X, and consumption of good X.")
print("\n--- Calculation Steps ---")
print(f"1. Change in nominal wage (w_hat): Based on the fixed price of traded good Y and fixed world return to capital, the equation is 0 = ({theta_ly:.3f})*w_hat + ({theta_ky:.3f})*0. This gives w_hat = {w_hat_pct:.4f}%.")
print(f"2. Change in price of X (Px_hat): The 2% tax on capital is passed through to the price. The equation is Px_hat = ({theta_lx:.3f})*({w_hat:.3f}) + ({theta_kx:.3f})*({tax_on_capital:.3f}). This gives Px_hat = {px_hat_pct:.4f}%.")
print(f"3. Change in income (I_hat): From the 3% increase in consumption of Y and an income elasticity of {eta_iy}, income change is {i_hat*100:.4f}%.")
print(f"4. Change in consumption of X (Cx_hat): Using the demand equation Cx_hat = ({epsilon_xx})*Px_hat + ({eta_ix})*I_hat, the calculation is Cx_hat = ({epsilon_xx})*({px_hat:.4f}) + ({eta_ix})*({i_hat:.4f}). This gives Cx_hat = {cx_hat_pct:.4f}%.")
print("\n--- Final Answer ---")
print("The percentage changes for (nominal wage, price of good X, consumption of good X) are:")
# Final comma-separated values
final_answer_string = f"{w_hat_pct},{px_hat_pct},{cx_hat_pct}"
print(final_answer_string)
print(f"<<<{final_answer_string}>>>")