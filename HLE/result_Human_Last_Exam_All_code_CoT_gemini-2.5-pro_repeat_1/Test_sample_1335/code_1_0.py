import math

# Step 0: Define the given parameters from the problem description.
# The "hat" notation (e.g., hat_w) will be used for percentage changes in decimal form.

# Tax on return to capital in sector X is 2%
tax_on_capital_X = 0.02

# The consumption of Y increases by 3%
hat_C_Y = 0.03

# Labor's share of cost in sector X and Y
theta_LX = 2/3
theta_LY = 1/3

# Capital's share of cost in sector X and Y
theta_KX = 1 - theta_LX  # This must be 1/3
theta_KY = 1 - theta_LY  # This must be 2/3

# Own-price elasticity of demand for good X and Y
eps_XX = -2.0
eps_YY = -1.0

# Income elasticity of demand for both goods
eta_IX = 1.0
eta_IY = 1.0

# Step 1: Calculate the percentage change in the nominal wage (hat_w).
# In a small open economy, the price of the traded good Y (P_Y) and the world return on capital (r) are fixed.
# This means their percentage changes, hat_P_Y and hat_r, are zero.
# The zero-profit condition in sector Y is: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r
# 0 = (1/3) * hat_w + (2/3) * 0
# This implies hat_w must be 0.
hat_w = 0.0
print("1. Calculation for the percentage change in nominal wage (w):")
print(f"The formula is: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r")
print(f"Plugging in the values: 0 = {theta_LY:.3f} * hat_w + {theta_KY:.3f} * 0")
print(f"This directly implies that hat_w = {hat_w}\n")

# Step 2: Calculate the percentage change in the price of good X (hat_P_X).
# The tax increases the cost of capital for producers in sector X by 2%.
hat_r_X = tax_on_capital_X
# The zero-profit condition in sector X is: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
print("2. Calculation for the percentage change in the price of good X (P_X):")
print(f"The formula is: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X")
print(f"Plugging in the values: hat_P_X = {theta_LX:.3f} * {hat_w} + {theta_KX:.3f} * {hat_r_X}")
print(f"The resulting change in the price of X is hat_P_X = {hat_P_X:.4f}\n")

# Step 3: Calculate the percentage change in nominal income (hat_I).
# We use the given information about the change in consumption of Y.
# The formula for the change in demand for Y is: hat_C_Y = eps_YY*hat_P_Y + eps_YX*hat_P_X + eta_IY*hat_I
# From the homogeneity property of demand functions: eps_YY + eps_YX + eta_IY = 0
eps_YX = -eps_YY - eta_IY
# Now we can solve for hat_I.
# hat_I = (hat_C_Y - eps_YY*hat_P_Y - eps_YX*hat_P_X) / eta_IY
# Since P_Y is fixed, hat_P_Y = 0.
hat_P_Y = 0
hat_I = (hat_C_Y - eps_YY * hat_P_Y - eps_YX * hat_P_X) / eta_IY
print("3. Calculation for the percentage change in nominal income (I):")
print(f"First, we find the cross-price elasticity eps_YX using the homogeneity property: eps_YY + eps_YX + eta_IY = 0")
print(f"Plugging in the values: {eps_YY} + eps_YX + {eta_IY} = 0  => eps_YX = {eps_YX}")
print(f"The formula for demand change is: hat_C_Y = eps_YY * hat_P_Y + eps_YX * hat_P_X + eta_IY * hat_I")
print(f"Plugging in the values: {hat_C_Y} = ({eps_YY}) * {hat_P_Y} + ({eps_YX}) * {hat_P_X:.4f} + ({eta_IY}) * hat_I")
print(f"Solving for hat_I gives: hat_I = {hat_I:.4f}\n")

# Step 4: Calculate the percentage change in the consumption of good X (hat_C_X).
# The formula for the change in demand for X is: hat_C_X = eps_XX*hat_P_X + eps_XY*hat_P_Y + eta_IX*hat_I
# From the homogeneity property of demand functions: eps_XX + eps_XY + eta_IX = 0
eps_XY = -eps_XX - eta_IX
# Now we can calculate hat_C_X.
hat_C_X = eps_XX * hat_P_X + eps_XY * hat_P_Y + eta_IX * hat_I
print("4. Calculation for the percentage change in consumption of good X (C_X):")
print(f"First, we find the cross-price elasticity eps_XY using the homogeneity property: eps_XX + eps_XY + eta_IX = 0")
print(f"Plugging in the values: {eps_XX} + eps_XY + {eta_IX} = 0  => eps_XY = {eps_XY}")
print(f"The formula for demand change is: hat_C_X = eps_XX * hat_P_X + eps_XY * hat_P_Y + eta_IX * hat_I")
print(f"Plugging in the values: hat_C_X = ({eps_XX}) * {hat_P_X:.4f} + ({eps_XY}) * {hat_P_Y} + ({eta_IX}) * {hat_I:.4f}")
print(f"The resulting change in the consumption of X is hat_C_X = {hat_C_X:.4f}\n")

# Step 5: Present the final answer as percentage changes, comma-separated.
dw_pct = hat_w * 100
dP_X_pct = hat_P_X * 100
dC_X_pct = hat_C_X * 100

print("--------------------------------------------------")
print("FINAL ANSWER")
print("The percentage changes for (1) nominal wage, (2) price of good X, and (3) consumption of good X are:")
print(f"{dw_pct},{dP_X_pct},{dC_X_pct}")
print("--------------------------------------------------")
<<<0.0,0.6666666666666666,1.6666666666666667>>>