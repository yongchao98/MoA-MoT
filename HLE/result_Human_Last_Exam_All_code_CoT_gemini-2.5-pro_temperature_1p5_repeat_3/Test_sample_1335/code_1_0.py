# Step 1: Define the given parameters in decimal form.
# Cost shares
theta_LX = 2/3  # Labor's share in X
theta_KX = 1 - theta_LX  # Capital's share in X
theta_LY = 1/3  # Labor's share in Y
theta_KY = 1 - theta_LY  # Capital's share in Y

# Changes in factor prices
r_hat_X = 0.02  # 2% tax on capital in sector X
r_hat_Y = 0.0   # No tax on capital in sector Y

# Changes in goods prices
P_y_hat = 0.0   # Y is freely traded (numeraire)

# Demand elasticities
epsilon_XX = -2 # Own-price elasticity for X
eta_X = 1       # Income elasticity for X

# Step 2: Calculate the percentage change in nominal wage (w_hat)
# From the zero-profit condition in sector Y: P_y_hat = theta_LY * w_hat + theta_KY * r_hat_Y
# 0 = (1/3) * w_hat + (2/3) * 0
w_hat = (-theta_KY * r_hat_Y - P_y_hat) / theta_LY
w_hat_pct = w_hat * 100

print("1. Calculation for percentage change in nominal wage (w_hat):")
print(f"   The equation is: P_y_hat = theta_LY * w_hat + theta_KY * r_y_hat")
print(f"   Substituting known values: {P_y_hat} = {theta_LY:.2f} * w_hat + {theta_KY:.2f} * {r_hat_Y}")
print(f"   Solving for w_hat gives: {w_hat:.4f}, or {w_hat_pct:.4f}%")
print("-" * 20)

# Step 3: Calculate the percentage change in the price of good X (Px_hat)
# From the zero-profit condition in sector X: Px_hat = theta_LX * w_hat + theta_KX * r_hat_X
Px_hat = theta_LX * w_hat + theta_KX * r_hat_X
Px_hat_pct = Px_hat * 100

print("2. Calculation for percentage change in price of good X (Px_hat):")
print(f"   The equation is: Px_hat = theta_LX * w_hat + theta_KX * r_hat_X")
print(f"   Substituting known values: Px_hat = {theta_LX:.2f} * {w_hat:.4f} + {theta_KX:.2f} * {r_hat_X}")
print(f"   Solving for Px_hat gives: {Px_hat:.4f}, or {Px_hat_pct:.4f}%")
print("-" * 20)

# Step 4: Calculate the percentage change in consumption of good X (Cx_hat)
# The key challenge is determining the change in income (I_hat). The provided static parameters are
# inconsistent. The simplest interpretation that allows a solution is that there is no change
# in the income of domestic consumers (e.g., capital is foreign-owned, tax not rebated). So, I_hat = 0.
I_hat = 0.0

# From the demand equation for X: Cx_hat = epsilon_XX * Px_hat + eta_X * I_hat
Cx_hat = epsilon_XX * Px_hat + eta_X * I_hat
Cx_hat_pct = Cx_hat * 100

print("3. Calculation for percentage change in consumption of good X (Cx_hat):")
print(f"   Assuming income change (I_hat) is 0.")
print(f"   The equation is: Cx_hat = epsilon_XX * Px_hat + eta_X * I_hat")
print(f"   Substituting known values: Cx_hat = {epsilon_XX} * {Px_hat:.4f} + {eta_X} * {I_hat}")
print(f"   Solving for Cx_hat gives: {Cx_hat:.4f}, or {Cx_hat_pct:.4f}%")
print("-" * 20)

# Final answer as comma-separated percentage values
print("Final answers as comma-separated percentage values (w_hat %, Px_hat %, Cx_hat %):")
print(f"{w_hat_pct},{Px_hat_pct},{Cx_hat_pct}")

# Final Answer Block
print("<<<" + f"{w_hat_pct},{Px_hat_pct},{Cx_hat_pct}" + ">>>")