# --- Model Parameters ---
# Factor cost shares in production
theta_LX = 2/3  # Labor's share in sector X
theta_KX = 1/3  # Capital's share in sector X
theta_LY = 1/3  # Labor's share in sector Y
theta_KY = 2/3  # Capital's share in sector Y

# Demand elasticities
eta_X = 1       # Income elasticity of demand for good X
eta_Y = 1       # Income elasticity of demand for good Y
epsilon_XX = -2 # Own-price elasticity of demand for good X

# --- Shocks and Known Values ---
# A 2% tax on capital in sector X leads to a 2% increase in its user cost
hat_r_X_prop = 0.02
# The consumption of good Y is given to increase by 3%
hat_C_Y_prop = 0.03
# Because the economy is small and open, the price of the traded good Y and the world rental rate on capital are fixed.
hat_P_Y_prop = 0
hat_r_Y_prop = 0

# --- Step-by-step Derivation ---
# The hat notation (e.g., hat_w) represents the proportional change (e.g., d_w / w).

# Step 1: Calculate the percentage change in nominal wage (w)
print("Step 1: Calculate the proportional change in nominal wage (hat_w).")
print("The zero-profit condition in sector Y is: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y")
print(f"Substituting known values: {hat_P_Y_prop} = (1/3) * hat_w + (2/3) * {hat_r_Y_prop}")
# This equation simplifies to 0 = (1/3) * hat_w, so hat_w must be 0.
hat_w_prop = 0
print(f"Result: hat_w = {hat_w_prop}\n")

# Step 2: Calculate the percentage change in the price of good X (P_X)
print("Step 2: Calculate the proportional change in the price of good X (hat_P_X).")
print("The zero-profit condition in sector X is: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X")
print(f"Substituting known values: hat_P_X = (2/3) * {hat_w_prop} + (1/3) * {hat_r_X_prop}")
hat_P_X_prop = theta_LX * hat_w_prop + theta_KX * hat_r_X_prop
print(f"Result: hat_P_X = {hat_P_X_prop}\n")

# Step 3: Calculate the percentage change in nominal income (I)
# The provided elasticities imply zero cross-price elasticities of demand.
print("Step 3: Calculate the proportional change in nominal income (hat_I).")
print("The demand function for good Y is: hat_C_Y = eta_Y * hat_I + ... (price terms)")
# Since hat_P_Y=0 and cross-price elasticity is 0, the equation simplifies.
print(f"Substituting known values into the simplified demand equation: {hat_C_Y_prop} = {eta_Y} * hat_I")
hat_I_prop = hat_C_Y_prop / eta_Y
print(f"Result: hat_I = {hat_I_prop}\n")

# Step 4: Calculate the percentage change in the consumption of good X (C_X)
print("Step 4: Calculate the proportional change in the consumption of good X (hat_C_X).")
print("The demand function for good X is: hat_C_X = eta_X * hat_I + epsilon_XX * hat_P_X")
print(f"Substituting calculated values: hat_C_X = {eta_X} * {hat_I_prop} + ({epsilon_XX}) * {hat_P_X_prop}")
hat_C_X_prop = eta_X * hat_I_prop + epsilon_XX * hat_P_X_prop
print(f"Result: hat_C_X = {hat_C_X_prop}\n")

# --- Final Answer ---
# Convert proportional changes to percentages by multiplying by 100
hat_w_perc = hat_w_prop * 100
hat_P_X_perc = hat_P_X_prop * 100
hat_C_X_perc = hat_C_X_prop * 100

print("---")
print("Final Answer: The percentage changes for (nominal wage, price of good X, consumption of good X) are:")
print(f"{hat_w_perc},{hat_P_X_perc},{hat_C_X_perc}")