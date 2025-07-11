import math

def combinations(n, k):
  """Calculates the number of combinations (n choose k)."""
  return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Number of points in the set P
num_points = 5
# The maximum number of points from P in any proper subcontinuum is k-1 = 2
max_points_in_subcontinuum = 2

# The largest number n is the number of ways to choose 2 points from 5.
n = combinations(num_points, max_points_in_subcontinuum)

# Print the equation and the result
print(f"The largest number n is C({num_points}, {max_points_in_subcontinuum}) = {n}")
print(f"Calculation: {num_points}! / ({max_points_in_subcontinuum}! * ({num_points}-{max_points_in_subcontinuum})!) = {math.factorial(num_points)} / ({math.factorial(max_points_in_subcontinuum)} * {math.factorial(num_points-max_points_in_subcontinuum)}) = {n}")