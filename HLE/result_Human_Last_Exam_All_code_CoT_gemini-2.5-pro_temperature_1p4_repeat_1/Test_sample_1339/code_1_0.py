import sys

def solve():
  """
  This function provides the solution to the theoretical group theory problem.
  
  (a) Existence and Uniqueness: The group described, hat(G), is the p-localization of G. 
  For a solvable group G, its p-localization exists and is unique up to isomorphism.
  
  (b) Maximum Derived Length: The group G has a subnormal series of length n with abelian factors,
  which implies its derived length dl(G) is at most n. The p-localization hat(G) also
  inherits a series of this type, with factors being the p-localization of the original abelian factors.
  These localized factors remain abelian. Therefore, hat(G) has a subnormal series of length n
  with abelian factors, which implies dl(hat(G)) <= n.
  This maximum can be achieved, for example, by considering groups of upper-triangular matrices
  over the integers, which can be constructed to have derived length n and satisfy the problem's
  hypotheses. For such groups, the derived length is preserved under p-localization.
  Thus, the maximum possible derived length is n.
  """
  
  # The answer to part (a) is Yes.
  part_a = "Yes"
  
  # The answer to part (b) is the expression "n".
  # It represents the length of the subnormal series given in the problem.
  part_b = "n"

  print(f"(a) {part_a}; (b) {part_b}")

solve()
