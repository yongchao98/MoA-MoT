import math

# Half-lives in days
T_Ba = 12.75
T_La = 1.678

# Measured activities in kBq/mL
A1 = 1.4
A2 = 2.1

# Time between measurements in days
dt = 14.0

# --- Step 1: Calculate decay constants (lambda) ---
# lambda = ln(2) / T_half-life
lambda_Ba = math.log(2) / T_Ba
lambda_La = math.log(2) / T_La

print(f"Decay constant for Ba-140 (lambda_Ba): {lambda_Ba:.5f} day^-1")
print(f"Decay constant for La-140 (lambda_La): {lambda_La:.5f} day^-1")
print("-" * 20)

# --- Step 2: Set up the ratio equation to solve for t1 ---
# The ratio of activities A2/A1 can be expressed as:
# A2/A1 = (exp(-lambda_Ba * (t1 + dt)) - exp(-lambda_La * (t1 + dt))) / (exp(-lambda_Ba * t1) - exp(-lambda_La * t1))
# This can be rearranged to solve for t1:
# t1 = ln(Z) / (lambda_La - lambda_Ba)
# where Z = (ratio_A - c2) / (ratio_A - c1)
# and c1 = exp(-lambda_Ba * dt), c2 = exp(-lambda_La * dt)

ratio_A = A2 / A1
print(f"Ratio of activities (A2/A1): {ratio_A:.3f}")

c1 = math.exp(-lambda_Ba * dt)
c2 = math.exp(-lambda_La * dt)

print(f"Decay factor for Ba-140 over {dt} days (c1): {c1:.5f}")
print(f"Decay factor for La-140 over {dt} days (c2): {c2:.5f}")

# --- Step 3: Solve for t1 ---
# t1 is the time between separation and the first measurement
numerator_Z = ratio_A - c2
denominator_Z = ratio_A - c1
Z = numerator_Z / denominator_Z

# The argument of the logarithm must be positive.
if Z <= 0:
    print("\nError: Mathematical inconsistency in the problem data. Cannot solve.")
else:
    t1 = math.log(Z) / (lambda_La - lambda_Ba)

    # --- Final Answer ---
    print("-" * 20)
    print(f"The calculated time between the chemical separation and the first analysis is {t1:.2f} days.")
    print("\nAssuming the analysis procedure started immediately after irradiation, this is the requested time.")
    final_equation = f"t = ln(({ratio_A:.3f} - {c2:.5f}) / ({ratio_A:.3f} - {c1:.5f})) / ({lambda_La:.5f} - {lambda_Ba:.5f}) = {t1:.2f} days"
    print("\nFinal calculation:")
    print(final_equation)

    # Reconstruct the equation with numbers to show the steps as requested by the user prompt
    print("\nThe time between sample irradiation and the first analysis is the solution `t` to the equation:")
    print(f"  {A2}/{A1} = (exp(-{lambda_Ba:.5f} * (t + {dt})) - exp(-{lambda_La:.5f} * (t + {dt}))) / (exp(-{lambda_Ba:.5f} * t) - exp(-{lambda_La:.5f} * t))")
    print(f"Solving this yields t = {t1:.2f} days.")
    
    # <<<1.03>>> is derived from the rounded output. 
    # The more precise result from the script is ~1.034. Let's provide the rounded answer.
    # The thought process resulted in 1.03.
    # final calculation: math.log(1.4493) / 0.3587 = 0.3711/0.3587 = 1.034
    
<<<1.03>>>