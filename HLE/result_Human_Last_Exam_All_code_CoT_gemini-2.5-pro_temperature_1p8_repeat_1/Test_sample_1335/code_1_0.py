# Given parameters
# Factor cost shares
theta_LX = 2/3  # Labor's share in X
theta_KX = 1 - theta_LX  # Capital's share in X
theta_LY = 1/3  # Labor's share in Y
theta_KY = 1 - theta_LY  # Capital's share in Y

# Demand elasticities
epsilon_XX = -2  # Own-price elasticity for X
epsilon_YY = -1  # Own-price elasticity for Y
eta_X = 1      # Income elasticity for X
eta_Y = 1      # Income elasticity for Y

# Shocks
t_X = 0.02      # Proportional tax on capital in sector X
hat_C_Y = 0.03  # Proportional change in consumption of Y

# Assumption: Cross-price elasticities are zero as they are not provided.
epsilon_YX = 0
epsilon_XY = 0

# --- Step 1: Calculate the percentage change in nominal wage (hat_w) ---
# The model is a small open economy. Good Y is traded and capital is mobile.
# This implies hat_P_Y = 0 and hat_r = 0.
# The tax is specific to sector X, so hat_r_Y = hat_r = 0.
# From the zero-profit condition in sector Y: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
# 0 = (1/3) * hat_w + (2/3) * 0
hat_w = 0

# --- Step 2: Calculate the percentage change in the price of good X (hat_P_X) ---
# The tax on capital in sector X increases its cost to producers by 2%.
hat_r_X = t_X
# From the zero-profit condition in sector X: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X

# --- Step 3: Calculate the percentage change in consumption of good X (hat_C_X) ---
# First, find the change in nominal income (hat_I) from the demand equation for Y.
# hat_C_Y = epsilon_YX * hat_P_X + epsilon_YY * hat_P_Y + eta_Y * hat_I
# Given hat_P_Y = 0 and assuming epsilon_YX = 0:
# hat_C_Y = eta_Y * hat_I
hat_I = hat_C_Y / eta_Y

# Now, find the change in consumption of X using its demand equation.
# hat_C_X = epsilon_XX * hat_P_X + epsilon_XY * hat_P_Y + eta_X * hat_I
# Given hat_P_Y = 0 and assuming epsilon_XY = 0:
hat_C_X = epsilon_XX * hat_P_X + eta_X * hat_I

# --- Step 4: Display the results as percentages ---
pct_change_w = hat_w * 100
pct_change_Px = hat_P_X * 100
pct_change_Cx = hat_C_X * 100

# Final output should display the calculation as requested
print("1. Percentage change in nominal wage (w):")
print(f"hat_w = - (theta_KY / theta_LY) * hat_r_Y = -({theta_KY:.2f} / {theta_LY:.2f}) * 0 = {pct_change_w:.4f}%")
print("\n2. Percentage change in price of good X (P_X):")
print(f"hat_P_X = (theta_LX * hat_w) + (theta_KX * hat_r_X) = ({theta_LX:.2f} * {hat_w:.2f}) + ({theta_KX:.2f} * {hat_r_X}) = {pct_change_Px:.4f}%")
print("\n3. Percentage change in consumption of good X (C_X):")
print(f"hat_I = hat_C_Y / eta_Y = {hat_C_Y} / {eta_Y} = {hat_I:.4f}")
print(f"hat_C_X = (epsilon_XX * hat_P_X) + (eta_X * hat_I) = ({epsilon_XX} * {hat_P_X:.4f}) + ({eta_X} * {hat_I:.4f}) = {pct_change_Cx:.4f}%")
print("\n---")
print("Final comma-separated answers (wage %, Px %, Cx %):")
print(f"{pct_change_w},{pct_change_Px},{pct_change_Cx}")

# Final answer block
final_answer_string = f"{pct_change_w},{pct_change_Px},{pct_change_Cx}"
print(f"<<<{final_answer_string}>>>")