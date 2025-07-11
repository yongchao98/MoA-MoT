import sys

# Step 0: Define all the given parameters as percentages
# Note: A hat (^) denotes a percentage change. For a variable V, hat_V = (dV/V)*100.

# Cost shares
theta_LX = 2/3  # Labor's share of cost in sector X
theta_KX = 1 - theta_LX  # Capital's share of cost in sector X
theta_LY = 1/3  # Labor's share of cost in sector Y
theta_KY = 1 - theta_LY  # Capital's share of cost in sector Y

# Demand Elasticities
eta_Px_X = -2  # Own-price elasticity of demand for X
eta_Py_Y = -1  # Own-price elasticity of demand for Y
eta_I_X = 1    # Income elasticity of demand for X
eta_I_Y = 1    # Income elasticity of demand for Y

# The shock
tax_rate = 2 # 2% tax
hat_rX_pct = tax_rate # The tax increases the cost of capital in sector X by 2%

# Known outcome
hat_CY_pct = 3 # Consumption of Y increases by 3%

# Since Y is a traded good and K is an internationally traded factor,
# their prices are fixed in the world market.
hat_PY_pct = 0
hat_r_pct = 0

# --- Calculations ---

# Step 1: Calculate the percentage change in nominal wage (hat_w)
# From the zero-profit condition in sector Y: hat_PY = theta_LY * hat_w + theta_KY * hat_r
# 0 = (1/3) * hat_w + (2/3) * 0
hat_w_pct = 0
print(f"Step 1: Calculate Percentage Change in Wage (w)")
print(f"Equation: 0 = ({theta_LY:.2f}) * hat_w + ({theta_KY:.2f}) * {hat_r_pct}")
print(f"Result: Percentage change in nominal wage = {hat_w_pct:.4f}%\n")


# Step 2: Calculate the percentage change in the price of good X (hat_PX)
# From the zero-profit condition in sector X: hat_PX = theta_LX * hat_w + theta_KX * hat_rX
# hat_PX = (2/3) * 0 + (1/3) * 2%
hat_PX_pct = theta_LX * hat_w_pct + theta_KX * hat_rX_pct
print(f"Step 2: Calculate Percentage Change in Price of X (Px)")
print(f"Equation: hat_Px = ({theta_LX:.2f}) * {hat_w_pct} + ({theta_KX:.2f}) * {hat_rX_pct}")
print(f"Result: Percentage change in price of X = {hat_PX_pct:.4f}%\n")


# Step 3: Calculate the percentage change in nominal income (hat_I)
# We use the demand function for good Y: hat_CY = eta_Px_Y * hat_PX + eta_Py_Y * hat_PY + eta_I_Y * hat_I
# First, find the cross-price elasticity eta_Px_Y using the homogeneity property of demand:
# eta_Px_Y + eta_Py_Y + eta_I_Y = 0
eta_Px_Y = -(eta_Py_Y + eta_I_Y)
# Now, solve for hat_I: hat_I = (hat_CY - eta_Px_Y * hat_PX - eta_Py_Y * hat_PY) / eta_I_Y
hat_I_pct = (hat_CY_pct - eta_Px_Y * hat_PX_pct - eta_Py_Y * hat_PY_pct) / eta_I_Y
print(f"Step 3: Calculate Percentage Change in Income (I)")
print(f"Derived Cross-Price Elasticity of Y (eta_Px_Y) from homogeneity: -({eta_Py_Y} + {eta_I_Y}) = {eta_Px_Y}")
print(f"Equation: {hat_CY_pct} = ({eta_Px_Y}) * ({hat_PX_pct:.2f}) + ({eta_Py_Y}) * {hat_PY_pct} + ({eta_I_Y}) * hat_I")
print(f"Result: Percentage change in nominal income = {hat_I_pct:.4f}%\n")


# Step 4: Calculate the percentage change in consumption of good X (hat_CX)
# We use the demand function for good X: hat_CX = eta_Px_X * hat_PX + eta_Py_X * hat_PY + eta_I_X * hat_I
# First, find the cross-price elasticity eta_Py_X using the homogeneity property of demand:
# eta_Px_X + eta_Py_X + eta_I_X = 0
eta_Py_X = -(eta_Px_X + eta_I_X)
# Now, calculate hat_CX:
hat_CX_pct = eta_Px_X * hat_PX_pct + eta_Py_X * hat_PY_pct + eta_I_X * hat_I_pct
print(f"Step 4: Calculate Percentage Change in Consumption of X (Cx)")
print(f"Derived Cross-Price Elasticity of X (eta_Py_X) from homogeneity: -({eta_Px_X} + {eta_I_X}) = {eta_Py_X}")
print(f"Equation: hat_Cx = ({eta_Px_X}) * ({hat_PX_pct:.2f}) + ({eta_Py_X}) * {hat_PY_pct} + ({eta_I_X}) * {hat_I_pct:.2f}")
print(f"Result: Percentage change in consumption of X = {hat_CX_pct:.4f}%\n")

# Final Answer
print("--- Final Answer ---")
print("The percentage changes for (nominal wage, price of good X, consumption of good X) are:")
print(f"{hat_w_pct},{hat_PX_pct},{hat_CX_pct}")
sys.stdout.flush()
# The final answer in the required format
# <<<0,0.6666666666666666,1.6666666666666667>>>