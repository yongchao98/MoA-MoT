import math

# --- Given Data ---
# Activity at first measurement (t=0) in kBq/mL
A1 = 1.4
# Activity at second measurement (t=14 days) in kBq/mL
A2 = 2.1
# Time between measurements in days
t_interval = 14
# Half-life of Yttrium-90 in hours
T_half_Y_hours = 64.1

# --- Step 1: Calculate the decay constant of Y-90 ---
# Convert half-life from hours to days
T_half_Y_days = T_half_Y_hours / 24
# Calculate decay constant (lambda) in units of days^-1
lambda_Y = math.log(2) / T_half_Y_days

# --- Step 2: Solve for the parent Sr-90 activity (A_Sr) ---
# The activity of Y-90 after 14 days is given by:
# A2 = A1 * exp(-lambda_Y * t_interval) + A_Sr * (1 - exp(-lambda_Y * t_interval))
# We rearrange to solve for A_Sr:
exp_term_14d = math.exp(-lambda_Y * t_interval)
A_Sr = (A2 - A1 * exp_term_14d) / (1 - exp_term_14d)

# --- Step 3: Solve for the time (T) from separation to the first measurement ---
# The initial measured activity A1 is due to the ingrowth of Y-90 over a time T
# A1 = A_Sr * (1 - exp(-lambda_Y * T))
# Rearranging for T: T = -ln(1 - A1/A_Sr) / lambda_Y
term_in_log = 1 - (A1 / A_Sr)
T_days = -math.log(term_in_log) / lambda_Y

# --- Final Output ---
print("This calculation determines the time elapsed between the chemical separation of Strontium-90 and the first analysis.")
print(f"The activity of the Sr-90 parent is calculated to be {A_Sr:.2f} kBq/mL.")
print("The time 'T' is found by solving the ingrowth equation:")
print(f"T = -ln(1 - A_Y(0) / A_Sr) / lambda_Y")
print(f"T = -ln(1 - {A1:.2f} / {A_Sr:.2f}) / {lambda_Y:.4f}")

# Final Result
print("\n--- Result ---")
print(f"The approximate time between separation and the first analysis is: {T_days:.2f} days.")
print(f"<<<{T_days:.2f}>>>")