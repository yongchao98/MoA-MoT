import sys

def f(n, m):
  """
  Determines if the first player has a winning position with a probability
  strictly more than 50% in 2D-Generalized NIM on a random n x m matrix.

  Args:
    n: The number of rows in the matrix.
    m: The number of columns in the matrix.

  Returns:
    1 if the first player's winning probability is > 50%, 0 otherwise.
  """
  # The logic is based on the analysis that only for the 1x1 case,
  # the probability of winning is not strictly greater than 50% (it is exactly 50%).
  # For all other board sizes, the number of P-positions (losing positions)
  # is significantly less than half the total number of positions.
  if n == 1 and m == 1:
    return 0
  else:
    return 1

def main():
  """
  Main function to demonstrate the usage of f(n, m).
  You can change the values of n and m here.
  """
  # Example values for n and m
  n = 2
  m = 2

  # The final result is determined by the logical equation: not (n == 1 and m == 1)
  # The numbers involved in this equation are n, m, and 1.
  # For n=2, m=2: not (2 == 1 and 2 == 1) -> not (False) -> True -> result is 1
  # For n=1, m=1: not (1 == 1 and 1 == 1) -> not (True) -> False -> result is 0
  
  result = f(n, m)
  
  print(f"For a game on a {n}x{m} matrix:")
  print(f"The value of f({n}, {m}) is: {result}")

if __name__ == "__main__":
  main()