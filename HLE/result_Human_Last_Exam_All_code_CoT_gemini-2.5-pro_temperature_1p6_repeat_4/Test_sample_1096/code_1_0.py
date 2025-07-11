import math

# Set the condition number M of the function.
# We assume m=1, so M is the condition number kappa.
M = 10.0

print(f"Calculating optimal gamma_1 and gamma_2 for a condition number M = {M}\n")

# The formulas for the optimal gamma values are:
# gamma_1,2 = ( 4*(M+1) +/- 2*sqrt(2)*(M-1) ) / (M**2 + 6*M + 1)
# We will calculate the numerator and denominator parts.

# Denominator D
D = M**2 + 6*M + 1

# Numerator parts
N1 = 4 * (M + 1)
N2 = 2 * math.sqrt(2) * (M - 1)

# Calculate gamma_1 and gamma_2
gamma_1 = (N1 + N2) / D
gamma_2 = (N1 - N2) / D

# Print the numbers in the final equations
print("The calculation is based on the formula: gamma = (N1 +/- N2) / D")
print(f"Denominator D = M**2 + 6*M + 1 = {M**2} + {6*M} + {1} = {D}")
print(f"Numerator term N1 = 4*(M+1) = 4*({M+1}) = {N1}")
print(f"Numerator term N2 = 2*sqrt(2)*(M-1) = 2*{math.sqrt(2)}*({M-1}) = {N2}")
print("-" * 30)

# Print the final results
print(f"The optimal choice for the pair (gamma_1, gamma_2) is:")
print(f"gamma_1 = ({N1} + {N2}) / {D} = {gamma_1}")
print(f"gamma_2 = ({N1} - {N2}) / {D} = {gamma_2}")
