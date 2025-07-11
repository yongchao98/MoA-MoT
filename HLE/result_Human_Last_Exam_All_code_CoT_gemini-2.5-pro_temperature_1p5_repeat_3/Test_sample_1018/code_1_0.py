import math

def simplicial_volume_surface(g):
  """
  Calculates the simplicial volume of a closed orientable surface of genus g.
  The formula is ||Σ_g|| = 4g - 4 for g >= 2, and 0 otherwise.
  """
  if g >= 2:
    return 4 * g - 4
  else:
    return 0

def combinations(n, k):
  """
  Calculates the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
  """
  if k < 0 or k > n:
    return 0
  return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Given genera and dimensions
g1 = 31
g2 = 17
n1 = 2
n2 = 2

# --- Step 1: Calculate simplicial volume of the first surface ---
sv1 = simplicial_volume_surface(g1)
print(f"The first surface is Σ_{g1} of genus g₁ = {g1}.")
print(f"Its simplicial volume is ||Σ_{g1}|| = 4 * {g1} - 4 = {sv1}.")
print("-" * 30)

# --- Step 2: Calculate simplicial volume of the second surface ---
sv2 = simplicial_volume_surface(g2)
print(f"The second surface is Σ_{g2} of genus g₂ = {g2}.")
print(f"Its simplicial volume is ||Σ_{g2}|| = 4 * {g2} - 4 = {sv2}.")
print("-" * 30)

# --- Step 3: Calculate the coefficient for the product formula ---
# The dimensions are n₁ = 2 and n₂ = 2.
n_total = n1 + n2
coeff = combinations(n_total, n1)
print(f"The product manifold is Σ_{g1} x Σ_{g2} of dimension {n1} + {n2} = {n_total}.")
print(f"The coefficient for the product formula is C({n1}+{n2}, {n1}) = C({n_total}, {n1}) = {coeff}.")
print("-" * 30)

# --- Step 4: Calculate the final simplicial volume of the product ---
total_sv = coeff * sv1 * sv2
print("The simplicial volume of the product is calculated as:")
print(f"||Σ_{g1} x Σ_{g2}|| = C({n_total}, {n1}) * ||Σ_{g1}|| * ||Σ_{g2}||")
print(f"||Σ_{31} x Σ_{17}|| = {coeff} * {sv1} * {sv2} = {total_sv}")
