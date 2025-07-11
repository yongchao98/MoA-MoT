import math

def calculate_max_zeros(n):
  """
  Calculates the maximal number of complex zeros for an N-channel problem.
  The formula is N * 2^(N-1).

  Args:
    n: The number of channels (N).

  Returns:
    The maximal number of zeros.
  """
  if not isinstance(n, int) or n <= 0:
    raise ValueError("N must be a positive integer.")
  
  power = n - 1
  result = n * (2**power)
  return result

if __name__ == "__main__":
  print("The maximal possible number of complex zeros is given by the formula Z = N * 2^(N-1).")
  print("Here are the calculations for some example values of N:")
  print("-" * 60)
  
  for n_val in range(1, 6):
    max_zeros = calculate_max_zeros(n_val)
    # The final equation is N * 2^(N-1) = Z
    # We print each number involved in this equation.
    print(f"For N = {n_val}, the number of zeros is {n_val} * 2^({n_val} - 1) = {max_zeros}")
