def print_treewidth_bound():
  """
  This function prints the formula for the tight upper bound on the treewidth
  of the union of two graphs H and G sharing k vertices.
  """
  t_H = "t_H"  # Treewidth of graph H
  t_G = "t_G"  # Treewidth of graph G
  k = "k"      # Number of shared vertices
  
  # The tight upper bound for the treewidth of the combined graph F is
  # max(t_H, t_G, k-1).
  
  print("The tight upper bound on the treewidth of F is given by the formula:")
  print(f"t_F <= max({t_H}, {t_G}, {k} - 1)")

print_treewidth_bound()