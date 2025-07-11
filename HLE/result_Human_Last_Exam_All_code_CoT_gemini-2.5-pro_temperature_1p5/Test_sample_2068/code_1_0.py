# Problem parameters
n = 2025
A = 10**15
B = 10**20

# The expression for the final value, as derived from the analysis, is (A * B) / (n - 1).
# We substitute the given numbers into this formula.
numerator = A * B
denominator = n - 1
result = numerator / denominator

# The final equation is (10**15 * 10**20) / (2025 - 1).
# The code outputs the numbers used in this equation and the final result.
print(f"Numerator A * B = {A} * {B} = {numerator}")
print(f"Denominator n - 1 = {n} - 1 = {denominator}")
print("Final result:")
print(result)
