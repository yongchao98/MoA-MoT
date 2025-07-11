# Define the given parameters
theta_Kx = 1/3  # Capital's share of cost in sector X
tax_rate = 0.02  # Tax on return to capital in sector X
eta_xx = -2      # Own-price elasticity of demand for good X
eta_xI = 1       # Income elasticity of demand for good X
Cy_hat = 0.03    # Percentage change in consumption of Y

# --- Calculations ---

# 1. Percentage change in nominal wage (w_hat)
# As derived in the explanation, the wage is fixed by the world prices of good Y and capital.
w_hat = 0.0

# 2. Percentage change in the price of good X (Px_hat)
# Px_hat = theta_Lx * w_hat + theta_Kx * r_x_hat
# With w_hat = 0 and r_x_hat = tax_rate
Px_hat = theta_Kx * tax_rate

# 3. Percentage change in real income (I_real_hat)
# From the demand for good Y, Cy_hat = eta_yx * Px_hat + I_real_hat.
# From homogeneity, eta_yx = 0.
# So, I_real_hat = Cy_hat.
I_real_hat = Cy_hat

# 4. Percentage change in consumption of good X (Cx_hat)
# Cx_hat = eta_xx * Px_hat + eta_xI * I_real_hat
Cx_hat = eta_xx * Px_hat + eta_xI * I_real_hat

# Convert to percentage points for the final answer
w_hat_pct = w_hat * 100
Px_hat_pct = Px_hat * 100
Cx_hat_pct = Cx_hat * 100

# The problem asks for the percentage change on nominal wage, price of good X, and the consumption of good X.
# We present the final answer as three comma-separated values.
print(f"{w_hat_pct},{Px_hat_pct},{Cx_hat_pct}")
