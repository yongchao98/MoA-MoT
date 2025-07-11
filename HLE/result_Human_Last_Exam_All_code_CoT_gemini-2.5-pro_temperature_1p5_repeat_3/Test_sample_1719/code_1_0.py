def get_obstruction_groups(n, k):
  """
  Prints the list of obstruction groups.

  Args:
    n: The dimension of the homology sphere X is (n-1).
    k: The rank of the vector bundle is 2k.
  """
  print("The homotopy-theoretic obstructions are classified by elements in the following groups:")
  print(f"1. pi_1(SO({2*k}))")
  print(f"2. pi_{n+1}(SO({2*k})), where n = {n}")

# The user's problem has a homology (n-1)-sphere.
# The user needs to provide the value of n.
# As an example, let's assume n=4, meaning X is a homology 3-sphere.
# As an example, let's assume k=3, meaning the bundle is rank 6.
n_example = 4
k_example = 3
get_obstruction_groups(n_example, k_example)
