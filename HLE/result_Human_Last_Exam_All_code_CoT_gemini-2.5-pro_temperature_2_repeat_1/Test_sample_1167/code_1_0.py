# Final result calculation script
# The reasoning is based on principles from harmonic analysis as described above.

alpha_numerator = -3
alpha_denominator = 4
alpha_value = alpha_numerator / alpha_denominator

print("The problem asks for the exponent alpha in the upper bound for the measure of the set X.")
print(f"Based on the analysis, the upper bound is found to be of the form C * log(N) * N^alpha.")
print("The dominant power-law term determines the value of alpha.")
print("\nThe final equation for alpha is:")
print(f"alpha = numerator / denominator")
print(f"where numerator = {alpha_numerator}")
print(f"and denominator = {alpha_denominator}")
print("\nThis gives the final value:")
print(f"alpha = {alpha_value}")
