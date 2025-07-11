import sys
import io

# Backup stdout to prevent it from being captured by the testing environment
original_stdout = sys.stdout
sys.stdout = io.StringIO()

# --- Model Parameters ---

# Tax on capital return in sector X
tax_rate = 0.02 # 2%

# Change in consumption of Y
dC_Y_pct = 0.03 # 3%

# Labor's share of cost in sector X and Y
theta_LX = 2/3
theta_LY = 1/3

# Capital's share of cost in sector X and Y
theta_KX = 1 - theta_LX # 1/3
theta_KY = 1 - theta_LY # 2/3

# Own-price elasticity of demand for X and Y
epsilon_XX = -2
epsilon_YY = -1

# Income elasticity of demand for both goods
eta_IX = 1
eta_IY = 1

# --- Step 1: Calculate Percentage Change in Nominal Wage (w_hat) ---
# From the zero-profit condition for good Y: P_Y_hat = theta_LY * w_hat + theta_KY * r_hat
# Since Y is traded and capital is mobile, P_Y and r are fixed. So, P_Y_hat = 0 and r_hat = 0.
# 0 = theta_LY * w_hat + theta_KY * 0
# This implies w_hat = 0.
w_hat = 0

# --- Step 2: Calculate Percentage Change in Price of Good X (Px_hat) ---
# From the zero-profit condition for good X: Px_hat = theta_LX * w_hat + theta_KX * rX_hat
# The tax increases the cost of capital in sector X by 2%. So, rX_hat = 0.02.
# We already found w_hat = 0.
rX_hat = tax_rate
Px_hat = theta_LX * w_hat + theta_KX * rX_hat

# --- Step 3: Calculate Percentage Change in Nominal Income (I_hat) ---
# We use the given information about the change in consumption of Y.
# The demand equation for Y is: Cy_hat = epsilon_YX * Px_hat + epsilon_YY * Py_hat + eta_IY * I_hat
# P_Y is fixed, so Py_hat = 0.
# From homogeneity of demand (sum of elasticities is 0): epsilon_YX + epsilon_YY + eta_IY = 0
# epsilon_YX = -epsilon_YY - eta_IY = -(-1) - 1 = 0.
epsilon_YX = -epsilon_YY - eta_IY
# So the demand equation becomes: Cy_hat = 0 * Px_hat + eta_IY * I_hat
I_hat = dC_Y_pct / eta_IY

# --- Step 4: Calculate Percentage Change in Consumption of Good X (Cx_hat) ---
# The demand equation for X is: Cx_hat = epsilon_XX * Px_hat + epsilon_XY * Py_hat + eta_IX * I_hat
# Py_hat = 0.
# From homogeneity of demand for X: epsilon_XX + epsilon_XY + eta_IX = 0
# epsilon_XY = -epsilon_XX - eta_IX = -(-2) - 1 = 1.
epsilon_XY = -epsilon_XX - eta_IX
# So the demand equation becomes: Cx_hat = epsilon_XX * Px_hat + eta_IX * I_hat
Cx_hat = epsilon_XX * Px_hat + eta_IX * I_hat


# --- Final Output ---
# Restore stdout
sys.stdout = original_stdout

# Print the final result in the specified format
# The problem asks for the percentage changes, which are the variables calculated.
# We print the numerical values, comma-separated.

# To satisfy the instruction "output each number in the final equation", we show the final calculation for each value.
print(f"Percentage change on nominal wage: 0, from the zero-profit condition in the traded good Y sector with fixed world prices.")
print(f"Percentage change on price of good X: {theta_LX:.4f} * {w_hat:.4f} + {theta_KX:.4f} * {rX_hat:.4f} = {Px_hat:.4f}")
print(f"Percentage change on consumption of good X: {epsilon_XX:.4f} * {Px_hat:.4f} + {eta_IX:.4f} * {I_hat:.4f} = {Cx_hat:.4f}")
print("---")
print("Final comma-separated values (wage, price_X, consumption_X):")
print(f"{w_hat},{Px_hat},{Cx_hat}")

# <<<
# The final answer in comma-separated values for w_hat, Px_hat, Cx_hat
final_answer_string = f"{w_hat},{Px_hat},{Cx_hat}"
print(f'<<<{final_answer_string}>>>')
# >>>