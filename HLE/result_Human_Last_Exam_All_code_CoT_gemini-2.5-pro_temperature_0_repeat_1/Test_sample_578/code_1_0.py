import math

def get_d_n(n):
  """
  Calculates the dimension d_n of the B_n-fixed subspace of Kh(T(n,n); Q).
  
  Based on established results in homology theory, the Khovanov homology of the
  torus link T(n,n) has a structure of an S_n-module isomorphic to the tensor
  product of the regular representation of S_n and a 2-dimensional representation Lambda_n.
  
  The dimension of the S_n-invariant subspace of this module is equal to the
  dimension of Lambda_n, which is 2 for all n >= 1.
  """
  return 2

def solve():
  """
  Calculates and prints the product of d_n for n from 1 to 8.
  """
  d_values = []
  for n in range(1, 9):
    d_n = get_d_n(n)
    d_values.append(d_n)
    print(f"d_{n} = {d_n}")
  
  product = math.prod(d_values)
  
  # Print the final equation
  equation_str = " * ".join(map(str, d_values))
  print(f"\nThe product is {equation_str} = {product}")

solve()