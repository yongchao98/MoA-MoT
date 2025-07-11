# Plan: Calculate the number of involutions in PSU(4, 997) using the derived formula.

# 1. Set the value of q.
q = 997

# 2. Define the formula for the number of involutions.
# The number of involutions is given by (term1 + term2) / 2.
# term1 corresponds to elements whose square is I (excluding the center).
# term2 corresponds to elements whose square is -I.

# 3. Calculate term1 = q^4 * (q^2 - q + 1) * (q^2 + 1)
term1 = q**4 * (q**2 - q + 1) * (q**2 + 1)

# 4. Calculate term2 = q^2 * (q^3 + 1)
term2 = q**2 * (q**3 + 1)

# 5. Calculate the total numerator.
numerator = term1 + term2

# 6. Calculate the final result by dividing by 2.
result = numerator // 2

# 7. Print the steps of the calculation as requested.
print("The number of involutions in PSU(4, q) for q=997 is calculated as follows:")
print("Formula: ( q^4 * (q^2 - q + 1) * (q^2 + 1) + q^2 * (q^3 + 1) ) / 2")
print(f"Let q = {q}")
print("\nStep 1: Calculate the first term of the numerator")
print(f"term1 = {q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1)")
print(f"term1 = {term1}")
print("\nStep 2: Calculate the second term of the numerator")
print(f"term2 = {q}^2 * ({q}^3 + 1)")
print(f"term2 = {term2}")
print("\nStep 3: Sum the terms to get the full numerator")
print(f"Numerator = term1 + term2 = {term1} + {term2}")
print(f"Numerator = {numerator}")
print("\nStep 4: Divide the numerator by 2 to get the final result")
print(f"Number of involutions = {numerator} / 2")
print(f"Result = {result}")
