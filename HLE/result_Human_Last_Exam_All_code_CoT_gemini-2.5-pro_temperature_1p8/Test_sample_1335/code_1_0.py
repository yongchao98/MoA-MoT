import sys

# Step 0: Define the given parameters
# Factor cost shares
theta_LX = 2/3  # Labor's share in X
theta_KX = 1/3  # Capital's share in X
theta_LY = 1/3  # Labor's share in Y
theta_KY = 2/3  # Capital's share in Y

# Elasticities of demand
e_XX = -2      # Own price elasticity of X
e_YY = -1      # Own price elasticity of Y
eta_IX = 1     # Income elasticity of X
eta_IY = 1     # Income elasticity of Y

# The economic shock
tax_rate = 0.02 # 2% tax on capital in sector X
dr_X = tax_rate

# Given outcome
dC_Y = 0.03 # 3% increase in consumption of Y

# Because Good Y is freely traded and capital is internationally mobile,
# we can treat their prices as fixed on the world market.
dP_Y = 0
dr_Y = 0

# --- Calculations ---

# Step 1: Calculate the percentage change in nominal wage (w)
# From the zero-profit condition in sector Y: dP_Y = theta_LY * dw + theta_KY * dr_Y
# 0 = (1/3) * dw + (2/3) * 0
dw = 0
print("1. Percentage change in nominal wage (w):")
print(f"Equation: 0 = {theta_LY:.2f} * dw + {theta_KY:.2f} * 0")
print(f"Result: dw = {dw*100:.3f}%\n")

# Step 2: Calculate the percentage change in the price of Good X (P_X)
# From the zero-profit condition in sector X: dP_X = theta_LX * dw + theta_KX * dr_X
dP_X = theta_LX * dw + theta_KX * dr_X
print("2. Percentage change in price of good X (P_X):")
print(f"Equation: dP_X = {theta_LX:.2f} * {dw*100:.2f}% + {theta_KX:.2f} * {dr_X*100:.2f}%")
print(f"Result: dP_X = {dP_X*100:.3f}%\n")


# Step 3: Calculate the percentage change in consumption of good X (C_X)
# First, find the percentage change in nominal income (I) using the given dC_Y
# The change in consumption of Y is: dC_Y = e_YY*dP_Y + e_YX*dP_X + eta_IY*dI
# We need the cross-price elasticity e_YX. From demand theory, e_YX + e_YY + eta_IY = 0.
e_YX = -e_YY - eta_IY
# So, dC_Y = e_YY*dP_Y + (-e_YY - eta_IY)*dP_X + eta_IY*dI. Wait, this is wrong.
# The correct homogeneity condition is for ONE good's demand function: e.g., for C_Y, e_YX + e_YY_price + eta_IY = 0. No. It is sum over all prices: e_Yx+e_Yy+eta_Y = 0. Let's re-write.
# Homogeneity of demand for Y: e_YX + e_YY + eta_IY = 0.
e_YX = - e_YY - eta_IY

# So, dI is solved from: dC_Y = e_YY*dP_Y + e_YX*dP_X + eta_IY*dI
dI = (dC_Y - e_YY * dP_Y - e_YX * dP_X) / eta_IY

print("3. Percentage change in consumption of good X (C_X):")
print("3a. First, find percentage change in nominal income (I):")
print(f"Homogeneity implies e_YX = -e_YY - eta_IY = -({e_YY}) - ({eta_IY}) = {e_YX}")
print(f"Equation for income: {dC_Y*100:.2f}% = {e_YY}*{dP_Y*100:.2f}% + {e_YX}*{dP_X*100:.3f}% + {eta_IY}*dI")
print(f"Solving for dI: dI = {dI*100:.3f}%\n")

# Now, calculate the percentage change in C_X
# The change in consumption of X is: dC_X = e_XX*dP_X + e_XY*dP_Y + eta_IX*dI
# We don't need e_XY because it's multiplied by dP_Y = 0.
dC_X = e_XX * dP_X + eta_IX * dI
print("3b. Second, find percentage change in consumption of X:")
print(f"Equation: dC_X = {e_XX} * {dP_X*100:.3f}% + {eta_IX} * {dI*100:.3f}%")
print(f"Result: dC_X = {dC_X*100:.3f}%\n")


# Final Answer
# Return the results as 3 comma-separated values in percentage terms
dw_perc = dw * 100
dPx_perc = dP_X * 100
dCx_perc = dC_X * 100

print("Final answers (percentage changes for w, P_X, C_X):")
print(f"{dw_perc},{dPx_perc},{dCx_perc}")

# To conform to the final output format request.
final_answer_str = f"{dw_perc},{dPx_perc},{dCx_perc}"
sys.stdout.write(f'<<<{final_answer_str}>>>')