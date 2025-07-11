def solve_manifold_problem():
  """
  This function provides the solution to the specified topology problem.
  
  The problem asks for the lexicographically least tuple (a_1, b_1, ..., a_l, b_l),
  with l minimal, such that M(a_i, b_i) is not "full", yet their connect-sum is "full".

  1. A manifold M(a,b) = M(a) x M(b) is not "full" if max(a,b) >= 2. This is based on
     an interpretation of the "full" property related to Stiefel-Whitney classes, which
     must be assumed for the problem to have a non-trivial solution.

  2. The connect-sum of two identical non-full manifolds, M(a,b) # M(a,b), becomes "full"
     because the relevant characteristic classes add up to zero in Z_2 cohomology.
     This means the minimal number of manifolds, l, is 2.

  3. We need to find the lexicographically smallest pair (a,b) such that max(a,b) >= 2.
     Listing pairs in order: (0,0), (0,1), (0,2), ...
     The first pair that satisfies the condition is (0,2).

  4. Therefore, the desired list of manifolds is M(0,2) and M(0,2).
     The resulting tuple is (a_1, b_1, a_2, b_2).
  """
  
  a1 = 0
  b1 = 2
  a2 = 0
  b2 = 2
  
  # The problem asks for the tuple to be printed.
  # The final equation is the tuple itself.
  result_tuple = (a1, b1, a2, b2)
  
  # The output format requires printing each number in the final "equation" (the tuple).
  # We will format it as a flat tuple string as requested in the prompt examples.
  print(f"({','.join(map(str, result_tuple))})")

solve_manifold_problem()