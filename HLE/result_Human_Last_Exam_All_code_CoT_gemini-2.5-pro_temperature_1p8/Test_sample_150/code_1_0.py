import sys

def solve_bit_plus_plus():
  """
  Executes a Bit++ program and prints the final value of x.
  The initial value of x is 0.
  """
  try:
    # Read the number of statements from the first line of input.
    num_statements = int(sys.stdin.readline())
  except (IOError, ValueError):
    # Handle cases with no or invalid input for the number of statements.
    num_statements = 0
  
  x = 0
  
  # Loop through each statement.
  for _ in range(num_statements):
    # Read a statement from a new line and remove leading/trailing whitespace.
    statement = sys.stdin.readline().strip()
    
    # Check if the statement contains the increment operator.
    # The problem guarantees that each statement is valid and will contain
    # either '++' or '--'.
    if '+' in statement:
      x += 1
    else:
      x -= 1
      
  # Print the final value of x.
  print(x)

solve_bit_plus_plus()