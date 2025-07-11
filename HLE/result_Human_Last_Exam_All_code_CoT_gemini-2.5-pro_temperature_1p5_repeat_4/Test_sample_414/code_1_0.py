import math

# Half-lives of the parent (Barium-140) and daughter (Lanthanum-140)
T_parent_days = 12.75  # Half-life of Ba-140 in days
T_daughter_hours = 40.28 # Half-life of La-140 in hours

# Convert daughter half-life to days
T_daughter_days = T_daughter_hours / 24.0

# Calculate decay constants (lambda = ln(2) / T_half) in units of days^-1
lambda_parent = math.log(2) / T_parent_days
lambda_daughter = math.log(2) / T_daughter_days

# Calculate the time to reach maximum daughter activity (t_max)
# t_max = ln(lambda_daughter / lambda_parent) / (lambda_daughter - lambda_parent)
numerator = math.log(lambda_daughter / lambda_parent)
denominator = lambda_daughter - lambda_parent
t_max_days = numerator / denominator

# The equation used is t_max = (ln(λ_d / λ_p)) / (λ_d - λ_p)
# We print each part of the calculation.
print(f"Half-life of Ba-140 (T_p): {T_parent_days:.2f} days")
print(f"Half-life of La-140 (T_d): {T_daughter_days:.2f} days")
print(f"Decay constant of Ba-140 (λ_p): {lambda_parent:.4f} day^-1")
print(f"Decay constant of La-140 (λ_d): {lambda_daughter:.4f} day^-1")
print("\nCalculating the time to maximum La-140 activity (t_max):")
print(f"t_max = ln({lambda_daughter:.4f} / {lambda_parent:.4f}) / ({lambda_daughter:.4f} - {lambda_parent:.4f})")
print(f"t_max = ln({lambda_daughter / lambda_parent:.4f}) / {denominator:.4f}")
print(f"t_max = {numerator:.4f} / {denominator:.4f}")
print(f"t_max = {t_max_days:.2f} days")

print(f"\nThe approximate time is {t_max_days:.1f} days.")
