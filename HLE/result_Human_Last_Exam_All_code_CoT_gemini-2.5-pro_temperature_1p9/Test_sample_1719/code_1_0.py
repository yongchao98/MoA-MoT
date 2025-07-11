import sys

def get_obstruction_groups():
  """
  This function defines the homotopy-theoretic obstruction groups.
  The variables n and k are treated as mathematical symbols.
  n-1 is the dimension of the homology sphere X.
  2k is the rank of the vector bundle E.
  """
  
  # For 2k >= 3, pi_1(SO(2k)) is Z/2Z. For 2k=2, it is Z.
  # We will represent it symbolically.
  group1 = "pi_1(SO(2k))"
  
  # The second group's index depends on n, the dimension of X + 1.
  group2 = "pi_{n+1}(SO(2k))"
  
  return [group1, group2]

if __name__ == '__main__':
  # The question asks for a list of homology groups of X
  # and/or homotopy groups of SO(2k). Our derivation shows the obstructions
  # are given by two homotopy groups of SO(2k). The homological information
  # about X is used to determine the index of the second group.
  
  obstruction_groups = get_obstruction_groups()
  
  print("The list of obstruction groups is:")
  for group in obstruction_groups:
    print(f"- {group}")
