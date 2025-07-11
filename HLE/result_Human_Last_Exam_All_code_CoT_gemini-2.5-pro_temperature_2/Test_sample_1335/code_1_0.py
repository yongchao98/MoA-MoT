import math

# --- Problem Parameters ---
# Factor cost shares
theta_LX = 2/3
theta_KX = 1 - theta_LX
theta_LY = 1/3
theta_KY = 1 - theta_LY

# Factor allocation shares
lambda_LX = 3/4
lambda_LY = 1 - lambda_LX

# Elasticities (as given in the problem statement)
sigma_X = -2
sigma_Y = -1

# Shock
d_r_X = 0.02  # 2% tax on capital in sector X -> percentage change in r_X is 2% or 0.02

# Given outcome
d_C_Y = 0.03  # Consumption of Y increases by 3%

# --- Step-by-Step Solution ---

# 1. Calculate the percentage change in the nominal wage (w_hat)
# From the zero-profit condition in sector Y: d(P_Y) = theta_LY * d(w) + theta_KY * d(r_Y)
# Since Y is traded and K is traded internationally, d(P_Y) = 0 and d(r_Y) = 0.
# 0 = theta_LY * d(w) + theta_KY * 0
d_w = 0.0

# 2. Calculate the percentage change in the price of good X (P_X_hat)
# From the zero-profit condition in sector X: d(P_X) = theta_LX * d(w) + theta_KX * d(r_X)
d_P_X = theta_LX * d_w + theta_KX * d_r_X

# 3. Calculate the percentage change in the consumption of good X (C_X_hat)
# This involves finding the change in output X based on production-side shifts.

# Relationship between factor price changes and factor quantity ratios
d_KL_ratio_X = sigma_X * (d_w - d_r_X)  # This is d(K_X/L_X)
d_KL_ratio_Y = sigma_Y * (d_w - 0)      # This is d(K_Y/L_Y)

# Express output change (dY) in terms of output change (dX)
# The full derivation leads to: dY = -(lambda_LX / lambda_LY) * dX + ...
# Let's derive it systematically:
# From factor substitution: dK_X = dL_X + d_KL_ratio_X; dK_Y = dL_Y + d_KL_ratio_Y
# From output equations: dX = dL_X + theta_KX * d_KL_ratio_X; dY = dL_Y + theta_KY * d_KL_ratio_Y
# From labor market clearing: dL_Y = -(lambda_LX / lambda_LY) * dL_X = -3 * dL_X
# Substitute:
# dL_X = dX - theta_KX * d_KL_ratio_X
# dL_Y = -3 * (dX - theta_KX * d_KL_ratio_X)
# dY = -3 * (dX - theta_KX * d_KL_ratio_X) + theta_KY * d_KL_ratio_Y
dY = -3 * (1) - (-3 * (-theta_KX * d_KL_ratio_X)) + (theta_KY * d_KL_ratio_Y)
slope_YX = -3
intercept_YX = -3 * (-theta_KX * d_KL_ratio_X) + (theta_KY * d_KL_ratio_Y)
intercept_YX = 3 * theta_KX * sigma_X * (d_w - d_r_X) + theta_KY * sigma_Y * (d_w-0)
# Substituting values:
intercept_val = 3 * (1/3) * (-2) * (0 - 0.02) + (2/3)*(-1)*(0)
intercept_val = 1 * (-2) * (-0.02) = 0.04

# So, the production relationship is dY = -3*dX + 0.04
# Using the assumption that dY = d_C_Y = 0.03
# 0.03 = -3*dX + 0.04
dX = (0.03 - 0.04) / -3

# For a non-traded good, production change equals consumption change
d_C_X = dX

# Convert results to percentages for the final output
wage_change_pct = d_w * 100
price_X_change_pct = d_P_X * 100
consumption_X_change_pct = d_C_X * 100

print(f"Percentage change on nominal wage: {wage_change_pct:.4f}")
print(f"Percentage change on price of good X: {price_X_change_pct:.4f}")
print(f"Percentage change on consumption of good X: {consumption_X_change_pct:.4f}")
print("\nFinal result as 3 comma-separated values:")
print(f"{wage_change_pct},{price_X_change_pct},{consumption_X_change_pct}")
