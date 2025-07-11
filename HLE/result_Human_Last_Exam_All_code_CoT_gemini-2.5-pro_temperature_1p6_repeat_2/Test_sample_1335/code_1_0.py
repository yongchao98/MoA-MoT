import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer


# --- Define Given Parameters ---
# Tax on capital in sector X
tax_k_x_pct = 2.0

# Labor's share of cost in sector X
theta_lx = 2/3
# Capital's share of cost in sector X
theta_kx = 1/3

# Labor's share of cost in sector Y
theta_ly = 1/3
# Capital's share of cost in sector Y
theta_ky = 2/3

# Demand elasticities
eta_ix = 1.0  # Income elasticity for X
eta_iy = 1.0  # Income elasticity for Y
eta_xx = -2.0 # Own-price elasticity for X
eta_yy = -1.0 # Own-price elasticity for Y

# Known outcome from the tax
c_y_change_pct = 3.0

# --- Calculations ---
# Convert percentages to decimals for internal calculations
t_k_x = tax_k_x_pct / 100.0
c_y_change = c_y_change_pct / 100.0

# --- Step 1: Calculate Percentage Change in Nominal Wage (w_hat) ---
# In sector Y, Price = Unit Cost. In percentage changes (hats):
# P_y_hat = theta_ly * w_hat + theta_ky * r_y_hat
# Since Y is traded and capital in Y is untaxed, P_y_hat = 0 and r_y_hat = 0.
# 0 = (1/3) * w_hat + (2/3) * 0 => w_hat = 0
w_change = 0.0
w_change_pct = w_change * 100

print("1. Calculation of the percentage change in nominal wage (w_hat):")
print("   From the zero-profit condition in the traded goods sector Y:")
print(f"   P_y_hat = theta_ly * w_hat + theta_ky * r_y_hat")
print(f"   0% = {theta_ly:.3f} * w_hat + {theta_ky:.3f} * 0%")
print(f"   Therefore, the percentage change in wage is {w_change_pct:.4f}%\n")


# --- Step 2: Calculate Percentage Change in Price of Good X (p_x_hat) ---
# In sector X, Price = Unit Cost. In percentage changes:
# p_x_hat = theta_lx * w_hat + theta_kx * r_x_hat
# r_x_hat is the 2% tax rate.
r_x_change = t_k_x
p_x_change = theta_lx * w_change + theta_kx * r_x_change
p_x_change_pct = p_x_change * 100

print("2. Calculation of the percentage change in the price of good X (p_x_hat):")
print("   From the zero-profit condition in sector X:")
print(f"   p_x_hat = theta_lx * w_hat + theta_kx * r_x_hat")
print(f"   p_x_hat = {theta_lx:.3f} * {w_change_pct:.4f}% + {theta_kx:.3f} * {tax_k_x_pct:.1f}%")
print(f"   Therefore, the percentage change in the price of X is {p_x_change_pct:.4f}%\n")


# --- Step 3: Calculate Percentage Change in Consumption of Good X (c_x_hat) ---
# This is a multi-step process.

# 3a: Find change in income (i_hat) from Y's demand equation.
# c_y_hat = eta_yx * p_x_hat + eta_yy * p_y_hat + eta_iy * i_hat
# First, find eta_yx using homogeneity: eta_yx + eta_yy + eta_iy = 0
eta_yx = -(eta_yy + eta_iy)
# Now solve for i_hat. Note p_y_hat = 0.
# i_hat = (c_y_hat - eta_yx * p_x_hat) / eta_iy
i_change = (c_y_change - eta_yx * p_x_change) / eta_iy
i_change_pct = i_change * 100

print("3. Calculation of the percentage change in the consumption of good X (c_x_hat):")
print("   First, we find the percentage change in national income (i_hat).")
print("   From the demand equation for good Y:")
print("   c_y_hat = eta_yx * p_x_hat + eta_yy * p_y_hat + eta_iy * i_hat")
print(f"   Using homogeneity (eta_yx + eta_yy + eta_iy = 0), we find eta_yx = -({eta_yy:.1f} + {eta_iy:.1f}) = {eta_yx:.1f}")
print(f"   {c_y_change_pct:.1f}% = {eta_yx:.1f} * {p_x_change_pct:.4f}% + {eta_yy:.1f} * 0% + {eta_iy:.1f} * i_hat")
print(f"   Therefore, the percentage change in income is {i_change_pct:.4f}%\n")

# 3b: Find change in consumption of X (c_x_hat).
# c_x_hat = eta_xx * p_x_hat + eta_xy * p_y_hat + eta_ix * i_hat
# First, find eta_xy using homogeneity: eta_xx + eta_xy + eta_ix = 0
eta_xy = -(eta_xx + eta_ix)
# Now solve for c_x_hat. Note p_y_hat = 0.
c_x_change = eta_xx * p_x_change + eta_xy * 0 + eta_ix * i_change
c_x_change_pct = c_x_change * 100

print("   Now, we find the percentage change in consumption of X.")
print("   From the demand equation for good X:")
print("   c_x_hat = eta_xx * p_x_hat + eta_xy * p_y_hat + eta_ix * i_hat")
print(f"   Using homogeneity (eta_xx + eta_xy + eta_ix = 0), we find eta_xy = -({eta_xx:.1f} + {eta_ix:.1f}) = {eta_xy:.1f}")
print(f"   c_x_hat = {eta_xx:.1f} * {p_x_change_pct:.4f}% + {eta_xy:.1f} * 0% + {eta_ix:.1f} * {i_change_pct:.4f}%")
print(f"   Therefore, the percentage change in consumption of X is {c_x_change_pct:.4f}%")


# --- Final Output ---
# The problem requests the final answer as 3 comma-separated values.
# The calculated values are for nominal wage, price of good X, and consumption of good X.
final_answer_string = f"{w_change_pct},{p_x_change_pct},{c_x_change_pct}"

# Restore original stdout
sys.stdout = original_stdout
# Get the content of the buffer
detailed_explanation = output_buffer.getvalue()
# Close the buffer
output_buffer.close()

# Print the detailed explanation
print(detailed_explanation)

# Print the final answer in the required format
# The question asks for the comma-separated values for the percentage changes.
# Percentage change in nominal wage = 0.0%
# Percentage change in price of good X = 0.6667%
# Percentage change in consumption of good X = 1.6667%
#
print("\nFinal Answer")
print("The percentage changes for [nominal wage, price of good X, consumption of good X] are:")
print(f"<<<{w_change_pct:.4f},{p_x_change_pct:.4f},{c_x_change_pct:.4f}>>>")