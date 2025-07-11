def print_complexity():
  """
  Prints the computational time complexity of the IMAG(n, c) algorithm.
  """
  
  # The time complexity is derived based on the bit complexity model,
  # considering that the cost of arithmetic operations depends on the bit length of the numbers.
  # The formula O((log n)^2 + log n * log c) combines the complexity for
  # different regimes of n and c (e.g., n > c and n < c).
  complexity_formula = "O((log n)^2 + log n * log c)"
  
  print("The computational time complexity of the IMAG(n, c) algorithm is:")
  print(complexity_formula)

print_complexity()