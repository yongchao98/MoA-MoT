import sys

def solve():
  """
  Executes the Bit++ program and prints the final value of x.
  """
  try:
    # Read the number of statements from the first line of input.
    n = int(sys.stdin.readline())
  except (ValueError, IndexError):
    # If input is empty or invalid, assume 0 statements.
    n = 0

  # The initial value of the variable x is 0.
  x = 0

  # Process each statement.
  for _ in range(n):
    statement = sys.stdin.readline()
    # If a '+' character is in the statement, it's an increment operation.
    if '+' in statement:
      x += 1
    # Otherwise, it must be a decrement operation.
    else:
      x -= 1

  # Print the final value of x.
  print(x)

solve()