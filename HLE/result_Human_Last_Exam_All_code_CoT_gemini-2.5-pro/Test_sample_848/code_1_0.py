import math

# This script calculates the value of 10^4 * lim (F(N)/ln(N)) as N -> infinity.
# The limit L is derived from the structure of the solutions to the given Diophantine conditions.
# L = 2/ln(alpha_13) + 2/ln(alpha_5)
# where alpha_k is the growth rate of the solution sequence for a parameter k.

# The two families of solutions correspond to k=13 and k=5.
k1 = 13
k2 = 5

# The growth rate alpha_k is the largest root of r^2 - k*r + 1 = 0.
# alpha_k = (k + sqrt(k^2 - 4)) / 2

# Calculate the growth rate for the first family of solutions (k=13)
alpha1 = (k1 + math.sqrt(k1**2 - 4)) / 2

# Calculate the growth rate for the second family of solutions (k=5)
alpha2 = (k2 + math.sqrt(k2**2 - 4)) / 2

# Calculate the natural logarithms of the growth rates
log_alpha1 = math.log(alpha1)
log_alpha2 = math.log(alpha2)

# Calculate the two terms of the limit formula
term1 = 2 / log_alpha1
term2 = 2 / log_alpha2

# Calculate the limit L
L = term1 + term2

# Calculate the final value required by the problem
final_value = 10000 * L

# Output the components of the calculation as requested
print("The final value is derived from the formula: 10^4 * (2/ln(alpha1) + 2/ln(alpha2))")
print(f"\nFor the k=13 family:")
print(f"  alpha1 = (13 + sqrt(13^2 - 4))/2 = (13 + sqrt(165))/2 = {alpha1}")
print(f"  ln(alpha1) = {log_alpha1}")
print(f"  Contribution to L = 2 / {log_alpha1} = {term1}")
print(f"\nFor the k=5 family:")
print(f"  alpha2 = (5 + sqrt(5^2 - 4))/2 = (5 + sqrt(21))/2 = {alpha2}")
print(f"  ln(alpha2) = {log_alpha2}")
print(f"  Contribution to L = 2 / {log_alpha2} = {term2}")
print(f"\nTotal Limit L = {term1} + {term2} = {L}")
print(f"Final Value = 10000 * L = 10000 * {L} = {final_value}")
print(f"\nThe integer part of the final value is {int(final_value)}.")
