import sys

def main():
  """
  Prints the new rules added by the Knuth-Bendix completion algorithm,
  ordered by their LHS using the specified LPO.

  This solution is based on the assumption that the third rule in the
  problem description contains a typo and should be f(g(x), h(x)) -> h(x),
  as the original system is unorientable.
  """
  
  # The new rules are derived from resolving critical pairs.
  # The rules are ordered increasingly by their left-hand side (LHS)
  # using the lexicographic path ordering with precedence f < g < h.
  
  rule1 = "f(f(y,y),f(y,y)) -> f(y,y)"
  rule2 = "g(f(y,y)) -> g(y)"
  rule3 = "h(x) -> g(x)"
  
  # Combine the rules into a single string, separated by commas.
  final_answer = f"{rule1}, {rule2}, {rule3}"
  
  print(final_answer)

if __name__ == "__main__":
  main()