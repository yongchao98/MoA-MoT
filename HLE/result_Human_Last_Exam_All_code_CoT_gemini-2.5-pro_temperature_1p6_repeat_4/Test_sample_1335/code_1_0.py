import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Model Parameters ---
# Elasticity of substitution (Note: The negative sign in the prompt is non-standard.
# We assume the standard positive definition, so sigma_x = 2, sigma_y = 1.
# This parameter is not needed for the final answer but is included for completeness.)
sigma_x = 2
sigma_y = 1

# Income elasticity of demand
eta_x = 1
eta_y = 1

# Own price elasticity of demand
eps_xx = -2
eps_yy = -1

# Labor's share of cost
theta_lx = 2/3
theta_ly = 1/3

# Capital's share of cost
theta_kx = 1 - theta_lx # 1/3
theta_ky = 1 - theta_ly # 2/3

# Tax rate on capital in sector X
tax_on_capital_x = 0.02

# Given outcome: change in consumption of Y
c_y_hat = 0.03

# --- Step 1: Calculate Percentage Change in Nominal Wage (w_hat) ---
# In a small open economy where good Y and capital K are freely traded,
# the prices P_y and r are fixed on the world market.
# The zero-profit condition in the Y sector is: P_y = c(w, r).
# In percentage changes: p_y_hat = theta_ly * w_hat + theta_ky * r_hat
# Since p_y_hat = 0 and r_hat = 0, it must be that w_hat = 0.
w_hat = 0.0

print("1. Calculation of the percentage change in nominal wage (w_hat):")
print(f"From the zero-profit condition in sector Y (p_y_hat = theta_ly*w_hat + theta_ky*r_hat):")
print(f"0 = ({theta_ly:.2f})*w_hat + ({theta_ky:.2f})*0")
print(f"Therefore, w_hat = {w_hat}")
print("-" * 20)

# --- Step 2: Calculate Percentage Change in Price of Good X (p_x_hat) ---
# The zero-profit condition for sector X is P_x = c(w, r_x), where r_x is the
# rental rate paid by producers in X, which includes the tax.
# In percentage changes: p_x_hat = theta_lx * w_hat + theta_kx * r_x_hat
# We know w_hat = 0 and r_x_hat is the tax rate of 2% or 0.02.
p_x_hat = theta_lx * w_hat + theta_kx * tax_on_capital_x

print("2. Calculation of the percentage change in the price of good X (p_x_hat):")
print(f"From the zero-profit condition in sector X (p_x_hat = theta_lx*w_hat + theta_kx*r_x_hat):")
print(f"p_x_hat = ({theta_lx:.2f})*({w_hat}) + ({theta_kx:.2f})*({tax_on_capital_x})")
print(f"p_x_hat = {p_x_hat}")
print("-" * 20)

# --- Step 3: Calculate Percentage Change in Consumption of Good X (c_x_hat) ---
# First, we determine the change in real income (I_real_hat) from the given change in C_y.
# The demand function for Y is: c_y_hat = eps_yx*p_x_hat + eps_yy*p_y_hat + eta_y*I_real_hat
# From demand theory (homogeneity), eps_yx + eps_yy + eta_y = 0.
# So, eps_yx = -eps_yy - eta_y = -(-1) - 1 = 0.
# The demand equation simplifies: c_y_hat = eta_y * I_real_hat
I_real_hat = c_y_hat / eta_y

print("3. Calculation of the percentage change in consumption of good X (c_x_hat):")
print("First, find the change in real income (I_real_hat) from c_y_hat = eta_y*I_real_hat:")
print(f"{c_y_hat} = ({eta_y})*I_real_hat  => I_real_hat = {I_real_hat}")

# Now, we use the demand function for X: c_x_hat = eps_xx*p_x_hat + eta_x*I_real_hat
c_x_hat = eps_xx * p_x_hat + eta_x * I_real_hat
print("\nNext, find the change in consumption of X (c_x_hat = eps_xx*p_x_hat + eta_x*I_real_hat):")
print(f"c_x_hat = ({eps_xx})*({p_x_hat}) + ({eta_x})*({I_real_hat})")
print(f"c_x_hat = {c_x_hat}")
print("-" * 20)

# --- Final Answer ---
# The question asks for 3 comma-separated values: w_hat, p_x_hat, c_x_hat
final_answer_string = f"{w_hat},{p_x_hat},{c_x_hat}"
print("\nThe final answer as 3 comma-separated values is:")
print(final_answer_string)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)

final_answer_formatted_for_submission = f"<<<{w_hat},{p_x_hat},{c_x_hat}>>>"
print(final_answer_formatted_for_submission)