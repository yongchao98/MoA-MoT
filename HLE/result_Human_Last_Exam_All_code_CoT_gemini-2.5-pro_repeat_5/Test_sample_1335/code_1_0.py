import math

# 1. Define the given parameters as proportional changes
# Tax on capital in sector X
tax_rate = 0.02

# Percentage change in factor costs and prices
E_rX = tax_rate  # Change in rental cost for producers in Sector X
E_rY = 0         # Change in rental cost for producers in Sector Y (world rate is fixed)
E_PY = 0         # Change in price of Good Y (traded good)

# Cost shares
theta_LX = 2/3  # Labor's share in X
theta_KX = 1/3  # Capital's share in X
theta_LY = 1/3  # Labor's share in Y
theta_KY = 2/3  # Capital's share in Y

# Demand elasticities
eta_XX = -2     # Own-price elasticity of demand for X
eta_YY = -1     # Own-price elasticity of demand for Y
eta_XI = 1      # Income elasticity of demand for X
eta_YI = 1      # Income elasticity of demand for Y

# Known outcome
E_CY = 0.03     # Percentage change in consumption of Y

# 2. Calculate the percentage change in nominal wage (E_w)
# In sector Y, the zero-profit condition is E(P_Y) = theta_LY * E(w) + theta_KY * E(r_Y)
# 0 = (1/3) * E_w + (2/3) * 0
# This implies E_w = 0.
E_w = 0

# 3. Calculate the percentage change in the price of good X (E_PX)
# In sector X, the zero-profit condition is E(P_X) = theta_LX * E(w) + theta_KX * E(r_X)
E_PX = theta_LX * E_w + theta_KX * E_rX

# 4. Calculate the percentage change in income (E_I)
# First, determine the cross-price elasticity of demand for Y with respect to the price of X (eta_YX).
# The homogeneity condition for demand functions states: eta_YX + eta_YY + eta_YI = 0
eta_YX = -eta_YY - eta_YI

# Now, use the demand function for Y to solve for E_I:
# E_CY = eta_YY*E_PY + eta_YX*E_PX + eta_YI*E_I
# 0.03 = (-1)*0 + (0)*E_PX + (1)*E_I
# This simplifies to E_I = E_CY.
E_I = (E_CY - eta_YY * E_PY - eta_YX * E_PX) / eta_YI

# 5. Calculate the percentage change in consumption of good X (E_CX)
# The demand function for X is: E_CX = eta_XX*E_PX + eta_XY*E_PY + eta_XI*E_I
# Since E_PY = 0, the cross-price term eta_XY*E_PY drops out.
E_CX = eta_XX * E_PX + eta_XI * E_I

# 6. Convert results to percentages for the final answer
wage_change_pct = E_w * 100
priceX_change_pct = E_PX * 100
consumptionX_change_pct = E_CX * 100

# Print the final answer as three comma-separated values:
# % change in nominal wage, % change in price of X, % change in consumption of X
print(f"{wage_change_pct},{priceX_change_pct},{consumptionX_change_pct}")