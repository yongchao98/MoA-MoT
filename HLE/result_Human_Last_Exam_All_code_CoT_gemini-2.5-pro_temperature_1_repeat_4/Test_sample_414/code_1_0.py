import math

# Step 1: Define constants
T_Ba = 12.75  # Half-life of Ba-140 in days
T_La = 40.3 / 24.0  # Half-life of La-140 in days (40.3 hours)
A1 = 1.4  # Initial measured activity in kBq/mL
A2 = 2.1  # Measured activity after 14 days in kBq/mL
t_decay = 14.0  # Time between measurements in days

# Step 2: Calculate decay constants (lambda)
lambda_Ba = math.log(2) / T_Ba
lambda_La = math.log(2) / T_La

print("--- Step 1 & 2: Constants and Decay Constants ---")
print(f"Half-life of Ba-140 (T_Ba): {T_Ba:.2f} days")
print(f"Half-life of La-140 (T_La): {T_La:.3f} days")
print(f"Decay constant for Ba-140 (λ_Ba): {lambda_Ba:.5f} day^-1")
print(f"Decay constant for La-140 (λ_La): {lambda_La:.5f} day^-1\n")

# Step 3: Calculate the ratio of Cherenkov efficiencies (R_eff = ε_La / ε_Ba)
# A2/A1 = exp(-λ_Ba*t) + R_eff * [λ_La/(λ_La-λ_Ba)] * [exp(-λ_Ba*t) - exp(-λ_La*t)]
# We solve for R_eff.
activity_ratio = A2 / A1
exp_Ba_decay = math.exp(-lambda_Ba * t_decay)
exp_La_decay = math.exp(-lambda_La * t_decay)

# Daughter growth term (without R_eff)
# Note: A_La_grown / A_Ba_initial = [λ_La/(λ_La-λ_Ba)] * [exp(-λ_Ba*t) - exp(-λ_La*t)]
daughter_growth_factor = (lambda_La / (lambda_La - lambda_Ba)) * (exp_Ba_decay - exp_La_decay)

# Solve for R_eff
R_eff = (activity_ratio - exp_Ba_decay) / daughter_growth_factor

print("--- Step 3: Calculate Cherenkov Efficiency Ratio (R_eff) ---")
print("The ratio of measured activities A2/A1 is the sum of the decayed parent's signal and the grown-in daughter's signal.")
print(f"A2/A1 = {A2}/{A1} = {activity_ratio:.3f}")
print("Equation: A2/A1 = exp(-λ_Ba*t) + R_eff * [λ_La/(λ_La-λ_Ba)] * [exp(-λ_Ba*t) - exp(-λ_La*t)]")
print(f"Solving for R_eff: ({activity_ratio:.3f} - {exp_Ba_decay:.3f}) / {daughter_growth_factor:.3f}")
print(f"The ratio of efficiencies (R_eff = ε_La / ε_Ba) is: {R_eff:.3f}\n")


# Step 4: Calculate the cooling time (t_cool) for maximum Cherenkov activity
# The time for maximum measured activity is given by:
# t_cool = (1/(λ_La-λ_Ba)) * ln[ (R_eff * K * λ_La) / (λ_Ba * (1 + R_eff*K)) ]
# where K = λ_La / (λ_La - λ_Ba)
K = lambda_La / (lambda_La - lambda_Ba)
numerator = R_eff * K * lambda_La
denominator = lambda_Ba * (1 + R_eff * K)
t_cool = (1 / (lambda_La - lambda_Ba)) * math.log(numerator / denominator)

print("--- Step 4: Calculate Cooling Time (t_cool) ---")
print("We assume the separation was performed at the time (t_cool) when the total Cherenkov signal from the unseparated sample was at its maximum.")
print("The equation for this time is:")
print("t_cool = (1 / (λ_La - λ_Ba)) * ln((R_eff * K * λ_La) / (λ_Ba * (1 + R_eff * K)))")
print(f"where K = λ_La / (λ_La - λ_Ba) = {K:.3f}\n")

print("--- Step 5: Final Calculation ---")
print("Plugging in the numbers into the equation for t_cool:")
# The equation is broken down for clarity
part1 = 1 / (lambda_La - lambda_Ba)
# Numerator of ln
num_ln = R_eff * K * lambda_La
# Denominator of ln
den_ln = lambda_Ba * (1 + R_eff * K)

print(f"t_cool = (1 / ({lambda_La:.4f} - {lambda_Ba:.4f})) * ln(({R_eff:.3f} * {K:.3f} * {lambda_La:.4f}) / ({lambda_Ba:.4f} * (1 + {R_eff:.3f} * {K:.3f})))")
print(f"t_cool = ({part1:.4f}) * ln({num_ln:.4f} / {den_ln:.4f})")
print(f"t_cool = {part1:.4f} * ln({num_ln/den_ln:.4f})")
print(f"t_cool = {part1:.4f} * {math.log(num_ln/den_ln):.4f}")
print(f"\nThe approximate time between irradiation and the first analysis is: {t_cool:.2f} days.")

final_answer = round(t_cool, 1)
# The prompt asks for the final answer in a specific format.
# However, printing it directly as a string might be better here.
# For example: <<<4.6>>>
# Let's print the value and then add the tag.
print(f"\nFinal Answer: {final_answer}")