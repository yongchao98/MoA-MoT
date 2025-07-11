import sympy

# Define the symbol n for our calculation
n = sympy.Symbol('n', positive=True)

# The potential kernel a(z) is asymptotically proportional to log|z|.
# We define a function for the potential based on the magnitude of the point z.
# The proportionality constant cancels out in the final ratio, so we can omit it.
def potential(magnitude):
  """Calculates the potential, proportional to log(magnitude)."""
  return sympy.log(magnitude)

# The magnitude of the center of the target disk, (n,0), is n.
mag_target_center = n

# The target disk is A_n = {z : |z - (n,0)| < n^(1/3)}.
# The point in A_n closest to the origin is (n - n^(1/3), 0).
# Its magnitude is n - n^(1/3).
mag_closest_point = n - n**(1/3)

# The probability p_n is approximated by the ratio of the potential
# at the closest point to the potential at the center of the target.
# p_n ≈ a(z_min) / a((n,0))
a_z_min = potential(mag_closest_point)
a_z_n = potential(mag_target_center)
p_n_expression = a_z_min / a_z_n

print("--- Symbolic Calculation of the Limit ---")
print(f"The potential at the closest point, a(z_min), is proportional to: {a_z_min}")
print(f"The potential at the target's center, a((n,0)), is proportional to: {a_z_n}")
print(f"The expression for the probability p_n is approximately:")
print(f"p_n ≈ ({a_z_min}) / ({a_z_n})")

# Now, we calculate the limit of this expression as n approaches infinity.
limit_p_n = sympy.limit(p_n_expression, n, sympy.oo)

print("\nThe final equation is:")
# Using sympy.Eq to format the output as an equation
final_equation = sympy.Eq(sympy.Limit(p_n_expression, n, sympy.oo), limit_p_n)
print(final_equation)
