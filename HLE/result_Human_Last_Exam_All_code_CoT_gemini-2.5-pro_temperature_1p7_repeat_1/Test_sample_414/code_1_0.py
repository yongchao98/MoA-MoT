import math

# Half-lives in days
T_half_Ba = 12.75
T_half_La = 1.678

# Calculate decay constants (lambda = ln(2) / T_half)
lambda_Ba = math.log(2) / T_half_Ba
lambda_La = math.log(2) / T_half_La

# Activity ratio and time difference
activity_ratio = 2.1 / 1.4
time_diff = 14.0 # in days

# The equation to solve for T is derived from the ratio of activities:
# e^((lambda_La - lambda_Ba) * T) = [activity_ratio - e^(-time_diff * lambda_La)] / [activity_ratio - e^(-time_diff * lambda_Ba)]

# Calculate the components of the right-hand side (RHS) of the equation
exp_term_La = math.exp(-time_diff * lambda_La)
exp_term_Ba = math.exp(-time_diff * lambda_Ba)

numerator = activity_ratio - exp_term_La
denominator = activity_ratio - exp_term_Ba

RHS = numerator / denominator

# The exponent term on the left-hand side (LHS)
lambda_diff = lambda_La - lambda_Ba

# Solve for T
# (lambda_La - lambda_Ba) * T = ln(RHS)
# T = ln(RHS) / (lambda_La - lambda_Ba)
ln_RHS = math.log(RHS)
T = ln_RHS / lambda_diff

print("Calculation Steps:")
print(f"1. Half-life of 140-Ba (T_p): {T_half_Ba} days")
print(f"2. Half-life of 140-La (T_d): {T_half_La} days")
print(f"3. Decay constant of 140-Ba (λ_p): ln(2)/{T_half_Ba:.3f} = {lambda_Ba:.5f} days⁻¹")
print(f"4. Decay constant of 140-La (λ_d): ln(2)/{T_half_La:.3f} = {lambda_La:.5f} days⁻¹")
print(f"5. Time between measurements (Δt): {time_diff} days")
print(f"6. Ratio of measured activities (A₂/A₁): {2.1} / {1.4} = {activity_ratio}")
print("\nSolving the equation for T, the time from irradiation to the first measurement:")
print("T = ln([ (A₂/A₁) - e^(-λ_d*Δt) ] / [ (A₂/A₁) - e^(-λ_p*Δt) ]) / (λ_d - λ_p)")
print(f"T = ln([ {activity_ratio} - e^(-{lambda_La:.5f}*{time_diff}) ] / [ {activity_ratio} - e^(-{lambda_Ba:.5f}*{time_diff}) ]) / ({lambda_La:.5f} - {lambda_Ba:.5f})")
print(f"T = ln([ {activity_ratio} - {exp_term_La:.5f} ] / [ {activity_ratio} - {exp_term_Ba:.5f} ]) / {lambda_diff:.5f}")
print(f"T = ln([ {numerator:.5f} ] / [ {denominator:.5f} ]) / {lambda_diff:.5f}")
print(f"T = ln({RHS:.5f}) / {lambda_diff:.5f}")
print(f"T = {ln_RHS:.5f} / {lambda_diff:.5f}")

print(f"\nResult:")
print(f"The approximate time between irradiation and the first analysis is {T:.2f} days.")
