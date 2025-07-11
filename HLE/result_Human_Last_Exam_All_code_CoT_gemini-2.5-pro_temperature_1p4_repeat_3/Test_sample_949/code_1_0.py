def solve():
  """
  This function determines the truth value of the six statements.
  As explained in the analysis, all statements are false. The family of
  k-ary tree strings T(k,d) serves as a counterexample for all of them
  by choosing a sufficiently large branching factor k.
  The universal constant required by the Big-O notation cannot exist.
  """
  
  # 1. sum(log L) = O(sum(log D)) -> False
  # 2. sum(loglog L) = O(sum(loglog D)) -> False
  # 3. sum(log^5 L) = O(sum(log^5 D)) -> False
  # 4. sum(2^sqrt(log L)) = O(sum(2^sqrt(log D))) -> False
  # 5. sum(L^0.1) = O(sum(D^0.11)) -> False
  # 6. sum(L^0.25) = O(sum(D^0.5)) -> False
  
  answer = "FFFFFF"
  print(answer)

solve()