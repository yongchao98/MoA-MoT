def print_obstruction_groups():
  """
  This function prints the list of homotopy-theoretic obstruction groups
  for the given problem. The variables n and k are parameters from the problem
  statement, where X is a homology (n-1)-sphere and E is a rank 2k bundle.
  """
  
  # The first obstruction group comes from evaluating the paths at a point.
  group1 = "pi_1(SO(2k))"
  
  # The second obstruction group comes from the based homotopy problem.
  # It is a cohomology group of X with coefficients in a homotopy group of SO(2k).
  group2 = "H^{n-1}(X, pi_{n+1}(SO(2k)))"
  
  print("The list of homotopy-theoretic obstruction groups is:")
  print(f"1. {group1}")
  print(f"2. {group2}")

print_obstruction_groups()