import math

# Half-lives of the Ba-140 / La-140 pair
T_Ba = 12.75  # days
T_La_hours = 40.272 # hours
T_La = T_La_hours / 24.0  # Convert La half-life to days

# Calculate the decay constants (lambda) in units of 1/day
lambda_Ba = math.log(2) / T_Ba
lambda_La = math.log(2) / T_La

# The time to reach maximum daughter (La-140) activity from a pure parent (Ba-140) source
# is given by the formula: t_max = ln(lambda_La / lambda_Ba) / (lambda_La - lambda_Ba)
t_max = (math.log(lambda_La / lambda_Ba)) / (lambda_La - lambda_Ba)

# The problem implies this characteristic time is the answer.
# We will now print the full calculation step-by-step.
print(f"Identifying the nuclide pair as Ba-140/La-140.")
print(f"The half-life of the parent, Ba-140, is {T_Ba} days.")
print(f"The half-life of the daughter, La-140, is {T_La_hours} hours or {T_La:.3f} days.")
print(f"The decay constant for Ba-140 (λ_Ba) = ln(2) / {T_Ba:.2f} days = {lambda_Ba:.5f} days⁻¹.")
print(f"The decay constant for La-140 (λ_La) = ln(2) / {T_La:.3f} days = {lambda_La:.5f} days⁻¹.")
print("\nThe cooling time 'T' is assumed to be the time required for the daughter activity to reach its maximum.")
print("This time is calculated using the formula: T = ln(λ_La / λ_Ba) / (λ_La - λ_Ba)")
print(f"\nPlugging in the values:")
print(f"T = ln({lambda_La:.5f} / {lambda_Ba:.5f}) / ({lambda_La:.5f} - {lambda_Ba:.5f})")
print(f"T = ln({lambda_La/lambda_Ba:.3f}) / ({lambda_La - lambda_Ba:.5f})")
print(f"T = {math.log(lambda_La/lambda_Ba):.5f} / {lambda_La - lambda_Ba:.5f}")
print(f"T = {t_max:.3f} days")

# The final answer is the calculated time, rounded to one decimal place as it is an "approximate" time.
final_answer = round(t_max, 1)
print(f"\nThe approximate time between irradiation and analysis is {final_answer} days.")
# No need to return the value in the specified format, the user can see the final printed line.
# However, to strictly adhere to the final instruction format.
# <<<5.7>>>