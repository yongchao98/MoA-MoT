# Define the given parameters and initial values
# Note: Percentage changes from the problem (e.g., 2%) are converted to decimals (0.02) for calculations.

# Cost shares
theta_LX = 2.0 / 3.0  # Labor's share in sector X
theta_KX = 1.0 / 3.0  # Capital's share in sector X
theta_LY = 1.0 / 3.0  # Labor's share in sector Y
theta_KY = 2.0 / 3.0  # Capital's share in sector Y

# Tax and resulting changes
hat_rX = 0.02  # 2% tax on capital in X
hat_CY = 0.03  # 3% increase in consumption of Y

# Demand elasticities
epsilon_XX = -2.0  # Own-price elasticity for X
epsilon_YY = -1.0  # Own-price elasticity for Y
eta_X = 1.0       # Income elasticity for X
eta_Y = 1.0       # Income elasticity for Y

# Step 1: Calculate the percentage change in the nominal wage (hat_w)
# From the zero-profit condition for the traded good Y: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
# In a small open economy, hat_P_Y = 0 and hat_r_Y = 0.
# 0 = theta_LY * hat_w + 0  => hat_w = 0
hat_w = 0.0

# Step 2: Calculate the percentage change in the price of good X (hat_PX)
# From the zero-profit condition for X: hat_PX = theta_LX * hat_w + theta_KX * hat_rX
hat_PX = theta_LX * hat_w + theta_KX * hat_rX

# Step 3: Calculate the percentage change in nominal income (hat_I)
# We use the demand function for Y: hat_CY = epsilon_YX * hat_PX + epsilon_YY * hat_PY + eta_Y * hat_I
# First, find the cross-price elasticity epsilon_YX from the homogeneity condition:
# epsilon_YX + epsilon_YY + eta_Y = 0
epsilon_YX = -epsilon_YY - eta_Y
# Since hat_PY = 0, the demand equation simplifies.
# hat_CY = epsilon_YX * hat_PX + eta_Y * hat_I
# Rearranging to solve for hat_I:
hat_I = (hat_CY - epsilon_YX * hat_PX) / eta_Y

# Step 4: Calculate the percentage change in the consumption of good X (hat_CX)
# Using the demand function for X: hat_CX = epsilon_XX * hat_PX + eta_X * hat_I (since hat_PY = 0)
hat_CX = epsilon_XX * hat_PX + eta_X * hat_I

# Convert results to percentages for the final output
hat_w_perc = hat_w * 100
hat_PX_perc = hat_PX * 100
hat_CX_perc = hat_CX * 100

# Print the final comma-separated values: % change in wage, % change in price of X, % change in consumption of X
# The question requires printing each number in the final equation. 
# The variables below hold the final calculated values derived from the equations in the plan.
print(f"{hat_w_perc},{hat_PX_perc},{hat_CX_perc}")
>>> 0.0,0.6666666666666666,1.6666666666666667