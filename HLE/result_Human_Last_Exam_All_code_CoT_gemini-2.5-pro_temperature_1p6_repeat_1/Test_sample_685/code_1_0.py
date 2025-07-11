import sys

def solve_game_probability(n, m):
  """
  Determines if the first player has a winning position with a probability 
  strictly more than 50% for a random n x m binary matrix.

  Args:
    n: The number of rows in the matrix.
    m: The number of columns in the matrix.

  Returns:
    1 if the probability is > 50%, 0 otherwise.
  """
  # Based on the analysis, the probability of the first player winning is
  # strictly greater than 50% for all cases except for a 1x1 matrix.
  # For n=1, m=1: P(win) is exactly 0.5, so f(1,1)=0.
  # For n>1 or m>1: P(win) > 0.5, so f(n,m)=1.
  if n == 1 and m == 1:
    return 0
  else:
    return 1

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python your_script_name.py <n> <m>")
    sys.exit(1)
  
  try:
    n_val = int(sys.argv[1])
    m_val = int(sys.argv[2])
    if n_val <= 0 or m_val <= 0:
      raise ValueError("Dimensions n and m must be positive integers.")
  except ValueError as e:
    print(f"Error: {e}")
    sys.exit(1)

  result = solve_game_probability(n_val, m_val)
  print(result)
