import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


# --- 1. Define Given Parameters ---
# Tax on capital in sector X
r_x_hat = 0.02
# Change in consumption of Y
c_y_hat = 0.03

# Cost shares
theta_LX = 2/3
theta_KX = 1 - theta_LX
theta_LY = 1/3
theta_KY = 1 - theta_LY

# Labor allocation shares
lambda_LX = 3/4
lambda_LY = 1 - lambda_LX

# Demand elasticities
eta_XX = -2.0  # Own-price elasticity for X
eta_YY = -1.0  # Own-price elasticity for Y
eta_I_X = 1.0   # Income elasticity for X
eta_I_Y = 1.0   # Income elasticity for Y

# Fixed world prices (percentage changes are zero)
p_y_hat = 0.0
r_y_hat = 0.0

print("### Step-by-step Calculation ###\n")

# --- 2. Calculate Percentage Change in Wage (w_hat) ---
print("--- Calculating Percentage Change in Wage (w_hat) ---")
print("The zero-profit condition for the traded good Y is: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat")
print(f"Given P_Y_hat=0, theta_LY={theta_LY:.2f}, and r_Y_hat=0, we have:")
# Equation: 0 = theta_LY * w_hat + theta_KY * 0
w_hat = 0
print(f"0 = {theta_LY:.2f} * w_hat => w_hat = {w_hat:.4f}\n")


# --- 3. Calculate Percentage Change in Price of X (p_x_hat) ---
print("--- Calculating Percentage Change in Price of X (p_x_hat) ---")
print("The zero-profit condition for good X is: P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat")
print(f"Using w_hat={w_hat:.4f}, theta_LX={theta_LX:.2f}, theta_KX={theta_KX:.2f}, and r_X_hat={r_x_hat}:")
p_x_hat = theta_LX * w_hat + theta_KX * r_x_hat
print(f"P_X_hat = {theta_LX:.2f} * {w_hat:.4f} + {theta_KX:.2f} * {r_x_hat}")
print(f"P_X_hat = {p_x_hat:.4f}\n")


# --- 4. Determine Demand-Side Parameters ---
print("--- Determining Demand-Side Parameters ---")
# First, find expenditure shares s_X and s_Y
# We know L_X/L_Y = lambda_LX/lambda_LY and P_i*Q_i = w*L_i/theta_Li
# So, (P_X*Q_X)/(P_Y*Q_Y) = (L_X/L_Y) * (theta_LY/theta_LX)
ratio_L = lambda_LX / lambda_LY
ratio_theta = theta_LY / theta_LX
ratio_PQ = ratio_L * ratio_theta
s_X = ratio_PQ / (1 + ratio_PQ)
s_Y = 1 - s_X
print("1. Calculate expenditure shares (s_X, s_Y):")
print(f"   Value of output ratio (P_X*Q_X)/(P_Y*Q_Y) = (L_X/L_Y) * (theta_LY/theta_LX) = {ratio_L:.2f} * ({theta_LY:.2f}/{theta_LX:.2f}) = {ratio_PQ:.2f}")
print(f"   This implies expenditure share for X, s_X = {ratio_PQ:.2f} / (1 + {ratio_PQ:.2f}) = {s_X:.2f}")
print(f"   And expenditure share for Y, s_Y = {s_Y:.2f}\n")

# Second, find cross-price elasticities using Cournot aggregation
print("2. Calculate cross-price elasticities (eta_YX, eta_XY):")
# From s_X*eta_XX + s_Y*eta_YX = -s_X
eta_YX = (-s_X - s_X * eta_XX) / s_Y
print(f"   From Cournot aggregation (s_X*eta_XX + s_Y*eta_YX = -s_X):")
print(f"   eta_YX = (-{s_X:.2f} - {s_X:.2f} * {eta_XX}) / {s_Y:.2f} = {eta_YX:.2f}")

# From s_X*eta_XY + s_Y*eta_YY = -s_Y
eta_XY = (-s_Y - s_Y * eta_YY) / s_X
print(f"   From Cournot aggregation (s_X*eta_XY + s_Y*eta_YY = -s_Y):")
print(f"   eta_XY = (-{s_Y:.2f} - {s_Y:.2f} * {eta_YY}) / {s_X:.2f} = {eta_XY:.2f}\n")


# --- 5. Calculate Percentage Change in Income (i_hat) ---
print("--- Calculating Percentage Change in Income (i_hat) ---")
print("Using the demand function for Y: C_Y_hat = eta_YX * P_X_hat + eta_YY * P_Y_hat + eta_I_Y * i_hat")
# Equation: c_y_hat = eta_YX * p_x_hat + eta_YY * p_y_hat + eta_I_Y * i_hat
# We solve for i_hat
i_hat = (c_y_hat - eta_YX * p_x_hat - eta_YY * p_y_hat) / eta_I_Y
print(f"Given C_Y_hat={c_y_hat}, we rearrange to find i_hat:")
print(f"i_hat = ({c_y_hat} - {eta_YX:.2f} * {p_x_hat:.4f} - {eta_YY:.2f} * {p_y_hat}) / {eta_I_Y:.2f}")
print(f"i_hat = {i_hat:.4f}\n")


# --- 6. Calculate Percentage Change in Consumption of X (c_x_hat) ---
print("--- Calculating Percentage Change in Consumption of X (c_x_hat) ---")
print("Using the demand function for X: C_X_hat = eta_XX * P_X_hat + eta_XY * P_Y_hat + eta_I_X * i_hat")
c_x_hat = eta_XX * p_x_hat + eta_XY * p_y_hat + eta_I_X * i_hat
print(f"C_X_hat = {eta_XX:.2f} * {p_x_hat:.4f} + {eta_XY:.2f} * {p_y_hat} + {eta_I_X:.2f} * {i_hat:.4f}")
print(f"C_X_hat = {c_x_hat:.4f}\n")


# --- 7. Final Result ---
# Convert proportional changes to percentages
w_hat_pct = w_hat * 100
p_x_hat_pct = p_x_hat * 100
c_x_hat_pct = c_x_hat * 100

print("--- Final Results (as percentages) ---")
print(f"Percentage change in nominal wage: {w_hat_pct:.4f}%")
print(f"Percentage change in price of good X: {p_x_hat_pct:.4f}%")
print(f"Percentage change in consumption of good X: {c_x_hat_pct:.4f}%")
print("\nFinal answer as 3 comma-separated values:")
print(f"{w_hat_pct:.4f},{p_x_hat_pct:.4f},{c_x_hat_pct:.4f}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_string = captured_output.getvalue()
print(output_string)

final_answer_string = f"<<<{w_hat_pct:.4f},{p_x_hat_pct:.4f},{c_x_hat_pct:.4f}>>>"
print(final_answer_string)