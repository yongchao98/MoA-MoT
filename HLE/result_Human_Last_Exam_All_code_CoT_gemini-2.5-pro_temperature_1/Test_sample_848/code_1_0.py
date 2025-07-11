import math

# This script calculates the final answer based on the plan above.
# The final result is derived from the expression:
# 10^4 * (2/ln(alpha) + 2/ln(gamma))
# where alpha and gamma are the characteristic roots of the recurrence relations
# for the two families of solutions (k=13 and k=5).

# Constants derived from the problem analysis
k1 = 13
d1 = 165  # Discriminant for k=13: 13^2 - 4 = 169 - 4 = 165
k2 = 5
d2 = 21   # Discriminant for k=5: 5^2 - 4 = 25 - 4 = 21

term_multiplier = 10000
numerator = 2

# The characteristic roots alpha and gamma
alpha = (k1 + math.sqrt(d1)) / 2
gamma = (k2 + math.sqrt(d2)) / 2

# Calculate the limit L
limit_L = numerator / math.log(alpha) + numerator / math.log(gamma)

# Calculate the final value requested by the problem
final_value = term_multiplier * limit_L

print("The final result is calculated from the expression:")
print(f"{term_multiplier} * ({numerator}/ln(({k1}+sqrt({d1}))/2) + {numerator}/ln(({k2}+sqrt({d2}))/2))")
print(f"Value of alpha = ({k1} + sqrt({d1})) / 2 = {alpha}")
print(f"Value of gamma = ({k2} + sqrt({d2})) / 2 = {gamma}")
print(f"Value of ln(alpha) = {math.log(alpha)}")
print(f"Value of ln(gamma) = {math.log(gamma)}")
print(f"Contribution from k={k1} family: {numerator / math.log(alpha)}")
print(f"Contribution from k={k2} family: {numerator / math.log(gamma)}")
print(f"Value of the limit L = {limit_L}")
print(f"Final result = {term_multiplier} * {limit_L} = {final_value}")
print(f"The integer part is: {int(final_value)}")
