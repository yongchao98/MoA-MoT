import math

def compute_lcm(a, b):
  """Computes the least common multiple of two integers."""
  return abs(a * b) // math.gcd(a, b)

# The denominators derived from index theory for G2 bundles are 3 and 8.
# The order of the torsion subgroup is related to their least common multiple.
denominator1 = 3
denominator2 = 8

# Compute the result
result = compute_lcm(denominator1, denominator2)

# The final result is the order of the largest cyclic torsion subgroup.
# We print the equation for clarity.
print(f"lcm({denominator1}, {denominator2}) = {result}")