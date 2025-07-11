import math

# Step 1 & 2: Define constants and the model
# Half-lives
T_half_Ba = 12.75  # days
T_half_La = 40.27 / 24.0  # convert hours to days

# Measured activities (of La-140 only)
A_La_t1 = 1.4  # kBq/mL at the first measurement
A_La_t2 = 2.1  # kBq/mL at the second measurement
time_diff = 14.0  # days between measurements

# Step 3: Calculate decay constants (lambda)
lambda_Ba = math.log(2) / T_half_Ba
lambda_La = math.log(2) / T_half_La

# Step 4: Use measurement data to find the state of the sample at the first measurement (t1)
# The equation for the daughter (La) activity is:
# A_La(t2) = A_La(t1)*exp(-lambda_La*dt) + [lambda_La / (lambda_La - lambda_Ba)] * A_Ba(t1) * [exp(-lambda_Ba*dt) - exp(-lambda_La*dt)]
# We need to solve for A_Ba(t1).

# Pre-calculate terms for clarity
R = lambda_La / (lambda_La - lambda_Ba)
exp_La_decay = math.exp(-lambda_La * time_diff)
exp_Ba_decay = math.exp(-lambda_Ba * time_diff)

# Rearrange the equation to solve for A_Ba_t1
numerator = A_La_t2 - A_La_t1 * exp_La_decay
denominator = R * (exp_Ba_decay - exp_La_decay)
A_Ba_t1 = numerator / denominator

print(f"Based on the measurements, the activities at the time of the first analysis were:")
print(f"Activity of 140-La = {A_La_t1:.3f} kBq/mL")
print(f"Activity of 140-Ba = {A_Ba_t1:.3f} kBq/mL\n")

# Step 5: Calculate the time T since the end of irradiation (EOI)
# At EOI, we assume A_La was 0.
# The ratio of activities at time T is given by:
# A_La(T) / A_Ba(T) = [lambda_La / (lambda_La - lambda_Ba)] * [1 - exp(-(lambda_La - lambda_Ba)*T)]
# We can solve this for T.

# Rearrange to solve for T
# A_La_t1 / A_Ba_t1 = R * (1 - exp(-(lambda_La - lambda_Ba)*T))
# (A_La_t1 / A_Ba_t1) / R = 1 - exp(-(lambda_La - lambda_Ba)*T)
# exp(-(lambda_La - lambda_Ba)*T) = 1 - (A_La_t1 / A_Ba_t1) / R
# -(lambda_La - lambda_Ba)*T = ln(1 - (A_La_t1 / A_Ba_t1) / R)
# T = -ln(1 - (A_La_t1 / A_Ba_t1) / R) / (lambda_La - lambda_Ba)

ratio_A = A_La_t1 / A_Ba_t1
term_in_log = 1 - ratio_A / R
T = -math.log(term_in_log) / (lambda_La - lambda_Ba)

print("The final calculation to find the time (T) is based on the equation:")
print(f"A_La(T) / A_Ba(T) = [λ_La / (λ_La - λ_Ba)] * [1 - exp(-(λ_La - λ_Ba) * T)]")
print("Plugging in the calculated and known values:")
print(f"{A_La_t1:.3f} / {A_Ba_t1:.3f} = {R:.3f} * (1 - exp(-({lambda_La:.3f} - {lambda_Ba:.3f}) * T))")
print("\nSolving this equation for T gives:")
print(f"Approximate time between irradiation and first analysis: {T:.2f} days")
