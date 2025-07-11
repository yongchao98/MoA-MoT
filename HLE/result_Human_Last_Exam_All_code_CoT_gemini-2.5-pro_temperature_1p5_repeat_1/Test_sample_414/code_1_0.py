import math

# Half-lives in days
T_Ba = 12.75  # Half-life of Ba-140
T_La = 40.27 / 24.0  # Half-life of La-140 in days (40.27 hours)

# Decay constants (lambda = ln(2) / T_half)
lambda_Ba = math.log(2) / T_Ba
lambda_La = math.log(2) / T_La

# Activity measurements and their ratio
A1 = 1.4  # in kBq/mL
A2 = 2.1  # in kBq/mL
activity_ratio = A2 / A1

# Time between measurements in days
delta_t = 14

# The equation to solve for t is:
# exp((lambda_La - lambda_Ba) * t) = (activity_ratio - exp(-lambda_La * delta_t)) / (activity_ratio - exp(-lambda_Ba * delta_t))

# Calculate the terms on the right-hand side of the rearranged equation
exp_term_La = math.exp(-lambda_La * delta_t)
exp_term_Ba = math.exp(-lambda_Ba * delta_t)

numerator = activity_ratio - exp_term_La
denominator = activity_ratio - exp_term_Ba

# Calculate the right-hand side
RHS = numerator / denominator

# Now solve for t
# (lambda_La - lambda_Ba) * t = ln(RHS)
t = math.log(RHS) / (lambda_La - lambda_Ba)

print("This plan solves for 't', the time between chemical separation and the first analysis.")
print(f"Half-life of Ba-140 (T_Ba): {T_Ba:.2f} days")
print(f"Half-life of La-140 (T_La): {T_La:.3f} days")
print(f"Decay constant for Ba-140 (λ_Ba): {lambda_Ba:.4f} day^-1")
print(f"Decay constant for La-140 (λ_La): {lambda_La:.4f} day^-1")
print(f"Activity Ratio (A2/A1): {activity_ratio:.1f}")
print(f"Time between measurements (Δt): {delta_t} days")
print("\nThe calculation is based on the formula:")
print("t = ln((A2/A1 - exp(-λ_La * Δt)) / (A2/A1 - exp(-λ_Ba * Δt))) / (λ_La - λ_Ba)")
print("\nPlugging in the numbers:")
print(f"t = ln(({activity_ratio:.1f} - {exp_term_La:.4f}) / ({activity_ratio:.1f} - {exp_term_Ba:.4f})) / ({lambda_La:.4f} - {lambda_Ba:.4f})")
print(f"t = ln({numerator:.4f} / {denominator:.4f}) / {lambda_La - lambda_Ba:.4f}")
print(f"t = ln({RHS:.4f}) / {lambda_La - lambda_Ba:.4f}")
print(f"t = {math.log(RHS):.4f} / {lambda_La - lambda_Ba:.4f}")
print("\nFinal Result:")
print(f"The approximate time between separation and the first analysis is {t:.2f} days.")

# The final answer is t, rounded to one decimal place as per convention.
final_answer = round(t, 1)