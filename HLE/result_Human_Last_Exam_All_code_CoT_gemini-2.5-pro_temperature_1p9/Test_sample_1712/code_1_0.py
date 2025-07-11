import math

# The problem is to evaluate the sum S = sum_{n=2 to infinity} (-1)^n * H_n / n
# where H_n is the n-th harmonic number.

# Through mathematical derivation, the exact value of the sum is found to be:
# S = 1 - (pi^2 / 12) + (ln(2))^2 / 2

# We can write the final equation as S = T1 - T2 + T3
T1 = 1.0
T2 = math.pi**2 / 12
T3 = (math.log(2))**2 / 2

# Calculate the final value S
S = T1 - T2 + T3

# Output the explanation, the value of each term in the equation, and the final result.
print("The analytical result of the sum is given by the equation:")
print("S = 1 - (pi^2 / 12) + (ln(2))^2 / 2\n")

print("Let's calculate the numerical value of each term:")
print(f"Term 1 (T1): {T1}")
print(f"Term 2 (pi^2 / 12): {T2}")
print(f"Term 3 ((ln(2))^2 / 2): {T3}\n")

print("The final equation with numerical values is:")
print(f"S = {T1} - {T2} + {T3}\n")

print("The final numerical value of the sum is:")
print(f"S = {S}")
