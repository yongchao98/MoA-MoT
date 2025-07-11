# --- 1. Define given parameters ---
# Tax on return to capital in sector X (as a decimal)
tax_rate = 0.02

# Cost shares
# Labor's share of cost in sector X is 2/3
theta_LX = 2/3
# Capital's share of cost in sector X is 1/3
theta_KX = 1/3
# Labor's share of cost in sector Y is 1/3
theta_LY = 1/3
# Capital's share of cost in sector Y is 2/3
theta_KY = 2/3

# Demand elasticities
# Own price elasticity of demand for good X is -2
epsilon_XX = -2
# Own price elasticity of demand for good Y is -1
epsilon_YY = -1
# Income elasticity of demand for both goods is 1
eta_IX = 1
eta_IY = 1

# Observed outcome
# Consumption of Y increases by 3% (as a decimal)
hat_C_Y = 0.03

# --- 2. Calculate percentage change in nominal wage (hat_w) ---
# In a small open economy with mobile capital and a traded good (Y),
# prices P_Y and r are fixed from the world market. So, hat_P_Y = 0 and hat_r = 0.
# The zero-profit condition for good Y is: P_Y = a_LY * w + a_KY * r.
# In percentage changes: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r
# 0 = (1/3) * hat_w + (2/3) * 0
# This implies hat_w = 0.
hat_w = 0.0

# --- 3. Calculate percentage change on the price of good X (hat_P_X) ---
# The 2% tax increases the cost of capital for producers in sector X.
hat_r_X = tax_rate
# The zero-profit condition for good X is: P_X = a_LX * w + a_KX * r_X.
# In percentage changes: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X

# --- 4. Calculate percentage change in nominal income (hat_I) ---
# We use the given change in consumption of Y.
# Demand for Y: hat_C_Y = epsilon_YY * hat_P_Y + epsilon_YX * hat_P_X + eta_IY * hat_I.
# Since hat_P_Y = 0, this simplifies to: hat_C_Y = epsilon_YX * hat_P_X + eta_IY * hat_I.
# We find epsilon_YX from the homogeneity condition of demand: epsilon_YY + epsilon_YX + eta_IY = 0
# -1 + epsilon_YX + 1 = 0  => epsilon_YX = 0
epsilon_YX = 0.0
# Substituting into the demand for Y equation:
# hat_C_Y = (0) * hat_P_X + eta_IY * hat_I
# hat_I = hat_C_Y / eta_IY
hat_I = hat_C_Y / eta_IY

# --- 5. Calculate percentage change on the consumption of good X (hat_C_X) ---
# We use the demand function for good X.
# hat_C_X = epsilon_XX * hat_P_X + epsilon_XY * hat_P_Y + eta_IX * hat_I.
# Since hat_P_Y = 0, this simplifies to:
hat_C_X = epsilon_XX * hat_P_X + eta_IX * hat_I

# --- 6. Print the results as comma-separated percentage values ---
# The results are percentage changes for: nominal wage, price of good X, and consumption of good X.
w_percent_change = hat_w * 100
Px_percent_change = hat_P_X * 100
Cx_percent_change = hat_C_X * 100

print(f"{w_percent_change},{Px_percent_change},{Cx_percent_change}")