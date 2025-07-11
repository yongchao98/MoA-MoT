def solve():
  """
  This function prints the formula for the number of ways to tile the L-shaped region.
  The formula is derived by decomposing the shape and analyzing the tiling configurations
  at the central corner.
  """
  number_part = 2
  fib_part_1 = "F_{n-1}"
  fib_part_2 = "F_n"
  
  print(f"The number of ways is given by the formula:")
  print(f"{number_part} * {fib_part_1} * {fib_part_2}")

solve()