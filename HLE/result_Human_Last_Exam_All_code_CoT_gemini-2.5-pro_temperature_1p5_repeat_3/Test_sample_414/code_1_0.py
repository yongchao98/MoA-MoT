import math

# Step 1: Define the constants for the Ba-140/La-140 system.
# Half-lives in days
T_Ba = 12.752  # Parent (Ba-140)
T_La = 1.6781  # Daughter (La-140)

# Activity measurements from the problem
A1 = 1.4  # kBq/mL at time T
A2 = 2.1  # kBq/mL at time T + 14 days
activity_ratio = A2 / A1

# Step 2: Calculate the decay constants (lambda) in units of days^-1.
lambda_Ba = math.log(2) / T_Ba
lambda_La = math.log(2) / T_La

# Step 3: Set up and solve the equation for the time T.
# The governing equation based on the activity ratio is:
# ratio = (exp(-lambda_Ba*(T+14)) - exp(-lambda_La*(T+14))) / (exp(-lambda_Ba*T) - exp(-lambda_La*T))
# This can be algebraically rearranged to solve for T:
# T = (1 / (lambda_La - lambda_Ba)) * ln((ratio - exp(-lambda_La*14)) / (ratio - exp(-lambda_Ba*14)))

# We calculate each part of the rearranged equation.
C_La = math.exp(-lambda_La * 14)
C_Ba = math.exp(-lambda_Ba * 14)

numerator = activity_ratio - C_La
denominator = activity_ratio - C_Ba
ln_argument = numerator / denominator
ln_term = math.log(ln_argument)
lambda_diff = lambda_La - lambda_Ba

# Calculate the final time T
T = (1 / lambda_diff) * ln_term

# Step 4: Print the full equation with the calculated values and the final result.
print("The problem is solved by modeling the ingrowth of La-140 from its parent Ba-140.")
print("The activity of the daughter, A_d(t), is proportional to: exp(-λ_p*t) - exp(-λ_d*t)")
print("\nThe ratio of the two activity measurements is:")
print(f"Ratio = {A2} / {A1} = {activity_ratio}")
print("\nWe solve the following equation for T, the time between irradiation and the first measurement:")
print(f"{activity_ratio:.2f} = (exp(-{lambda_Ba:.5f} * (T + 14)) - exp(-{lambda_La:.5f} * (T + 14))) / (exp(-{lambda_Ba:.5f} * T) - exp(-{lambda_La:.5f} * T))")
print("\nThe rearranged equation used for calculation is:")
print("T = [1 / (λ_La - λ_Ba)] * ln[(Ratio - exp(-λ_La*14)) / (Ratio - exp(-λ_Ba*14))]")
print("\nPlugging in the numerical values:")
print(f"T = [1 / ({lambda_La:.5f} - {lambda_Ba:.5f})] * ln[({activity_ratio:.2f} - {C_La:.5f}) / ({activity_ratio:.2f} - {C_Ba:.5f})]")
print(f"T = [1 / {lambda_diff:.5f}] * ln[{numerator:.5f} / {denominator:.5f}]")
print(f"T = {1/lambda_diff:.5f} * {ln_term:.5f}")

print("\n-------------------------------------------------------------------------")
print(f"The approximate time between irradiation and the first analysis is: {T:.2f} days.")
print("-------------------------------------------------------------------------")