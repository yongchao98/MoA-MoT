import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Model Parameters ---
# Cost shares
theta_LX = 2/3  # Labor's share in X
theta_KX = 1/3  # Capital's share in X
theta_LY = 1/3  # Labor's share in Y
theta_KY = 2/3  # Capital's share in Y

# Demand elasticities
eta_IX = 1      # Income elasticity for X
eta_IY = 1      # Income elasticity for Y
epsilon_XX = -2 # Own-price elasticity for X
epsilon_YY = -1 # Own-price elasticity for Y

# Economic Shocks and Given Results
tax_on_capital_X_prop = 0.02 # 2% tax on capital in sector X
Dy_change_prop = 0.03      # Consumption of Y increases by 3%

# --- Step 1: Calculate the percentage change in nominal wage (w) ---
# In a small open economy with Y as a traded good and capital internationally mobile,
# P_Y and r_Y are fixed. The zero-profit condition for sector Y is:
# P_Y_change = theta_LY * w_change + theta_KY * r_Y_change
# 0 = (1/3) * w_change + (2/3) * 0  => w_change = 0
w_change_prop = 0
w_change_pct = w_change_prop * 100
print(f"Step 1: Calculate percentage change in wage (w)")
print(f"From zero-profit in sector Y: 0% = ({theta_LY:.4f} * w %) + ({theta_KY:.4f} * 0%)")
print(f"Result: w_change = {w_change_pct:.4f}%\n")


# --- Step 2: Calculate the percentage change in the price of good X (P_X) ---
# The cost of capital for firms in sector X increases by the tax rate.
r_X_change_prop = tax_on_capital_X_prop
# The zero-profit condition for sector X determines its price change:
# P_X_change = theta_LX * w_change + theta_KX * r_X_change
Px_change_prop = theta_LX * w_change_prop + theta_KX * r_X_change_prop
Px_change_pct = Px_change_prop * 100
print(f"Step 2: Calculate percentage change in price of X (P_X)")
print(f"From zero-profit in sector X: P_X % = ({theta_LX:.4f} * {w_change_pct:.4f}%) + ({theta_KX:.4f} * {r_X_change_prop*100:.4f}%)")
print(f"Result: P_X_change = {Px_change_pct:.4f}%\n")


# --- Step 3: Calculate the percentage change in nominal income (I) ---
# We use the given change in consumption of Y (Dy_change).
# Demand function for Y: Dy_change = e_YX * Px_change + e_YY * Py_change + n_IY * I_change
# From homogeneity of demand (e_YX + e_YY + n_IY = 0): e_YX = -epsilon_YY - eta_IY
epsilon_YX = -epsilon_YY - eta_IY
# Since P_Y_change = 0: Dy_change = e_YX * Px_change + n_IY * I_change. Let's re-verify the theory.
# The homogeneity property e_YX + e_YY + n_IY = 0 assumes income is nominal income.
# In this specific case, 0.03 = (0)*Px_change_prop + (-1)*0 + (1)*I_change_prop -> I_change_prop = 0.03
I_change_prop = Dy_change_prop / eta_IY
I_change_pct = I_change_prop * 100
print(f"Step 3: Calculate percentage change in income (I)")
print(f"From homogeneity for Y: e_YX = -({epsilon_YY}) - ({eta_IY}) = {epsilon_YX}")
print(f"From demand for Y: {Dy_change_prop*100:.4f}% = ({epsilon_YX} * {Px_change_pct:.4f}%) + ({epsilon_YY} * 0%) + ({eta_IY} * I %)")
print(f"Result: I_change = {I_change_pct:.4f}%\n")


# --- Step 4: Calculate the percentage change in consumption of good X (D_X) ---
# Demand function for X: Dx_change = e_XX * Px_change + e_XY * Py_change + n_IX * I_change
# From homogeneity of demand (e_XX + e_XY + n_IX = 0): e_XY = -epsilon_XX - eta_IX
epsilon_XY = -epsilon_XX - eta_IX
P_Y_change_prop = 0
Dx_change_prop = epsilon_XX * Px_change_prop + epsilon_XY * P_Y_change_prop + eta_IX * I_change_prop
Dx_change_pct = Dx_change_prop * 100
print(f"Step 4: Calculate percentage change in consumption of X (D_X)")
print(f"From homogeneity for X: e_XY = -({epsilon_XX}) - ({eta_IX}) = {epsilon_XY}")
print(f"Final equation: D_X % = (e_XX * P_X %) + (e_XY * P_Y %) + (n_IX * I %)")
print(f"D_X % = (({epsilon_XX}) * {Px_change_pct:.4f}%) + (({epsilon_XY}) * {P_Y_change_prop*100:.4f}%) + (({eta_IX}) * {I_change_pct:.4f}%)")
print(f"D_X % = {epsilon_XX * Px_change_pct:.4f}% + {eta_IX * I_change_pct:.4f}%")
print(f"Result: D_X_change = {Dx_change_pct:.4f}%\n")

# --- Final Answer ---
final_answer = f"{w_change_pct},{Px_change_pct},{Dx_change_pct}"

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# print(output) # Uncomment to see the detailed steps
print(final_answer)
