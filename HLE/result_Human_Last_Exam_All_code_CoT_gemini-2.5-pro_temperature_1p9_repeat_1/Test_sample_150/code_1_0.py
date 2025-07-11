import sys

def solve_bit_plus_plus():
  """
  Executes a Bit++ program and prints the final value of x.
  The initial value of x is 0.
  """
  try:
    # Read the number of statements from standard input.
    n = int(sys.stdin.readline())
  except (IOError, ValueError):
    # If input is empty or not a valid number, assume 0 statements.
    n = 0
  
  x = 0
  # Process each statement.
  for _ in range(n):
    statement = sys.stdin.readline()
    
    # A '+' in the statement means it's an increment operation.
    if '+' in statement:
      x += 1
    # Otherwise, it must be a decrement operation.
    else:
      x -= 1
      
  # Print the final value of x.
  print(x)

solve_bit_plus_plus()