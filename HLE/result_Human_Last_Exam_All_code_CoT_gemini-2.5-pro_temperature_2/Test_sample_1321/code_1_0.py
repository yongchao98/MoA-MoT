def a_1d(n):
  """
  Calculates the value of a(n) in the 1D case.
  The formula is a(n) = 2n - 2.
  """
  if n < 1:
    return 0
  return 2 * n - 2

# Calculate a(2)
n1 = 2
result1 = a_1d(n1)
print(f"For n = {n1}, the calculation is: a({n1}) = 2 * {n1} - 2 = {result1}")

# Calculate a(3)
n2 = 3
result2 = a_1d(n2)
print(f"For n = {n2}, the calculation is: a({n2}) = 2 * {n2} - 2 = {result2}")

# Calculate a(42)
n3 = 42
result3 = a_1d(n3)
print(f"For n = {n3}, the calculation is: a({n3}) = 2 * {n3} - 2 = {result3}")
