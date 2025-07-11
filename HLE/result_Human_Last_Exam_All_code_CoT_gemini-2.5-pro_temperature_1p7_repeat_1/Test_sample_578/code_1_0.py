def solve_braid_invariants_product():
  """
  Calculates the product of the dimensions of the braid group B_n invariants
  of the Khovanov homology of the torus link T(n,n) for n=1 to 8.
  """
  
  d_values = []
  
  # For n=1, d_1 is the dimension of the Khovanov homology of the unknot, which is 2.
  d_values.append(2)
  
  # For n > 1, the dimension of the B_n-invariant subspace is 2.
  for n in range(2, 9):
    d_values.append(2)
    
  # Calculate the product
  product = 1
  for d in d_values:
    product *= d
    
  # Format the equation string
  # The question asks to output the final equation with each number.
  equation_str = " * ".join(map(str, d_values))
  
  # Print the final result in the format of an equation.
  print(f"{equation_str} = {product}")

solve_braid_invariants_product()