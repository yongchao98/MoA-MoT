import sys

def solve_bit_plus_plus():
  """
  Solves the Bit++ problem by reading statements from standard input.
  The initial value of x is 0. It is increased for statements containing '++'
  and decreased for statements containing '--'.
  """
  try:
    # Read the number of statements.
    n = int(sys.stdin.readline())
    
    # Initialize the variable x to 0.
    x = 0
    
    # Process each statement.
    for _ in range(n):
      statement = sys.stdin.readline().strip()
      
      # The problem guarantees that each statement contains either '++' or '--'.
      # We can check for the presence of a '+' character to determine the operation.
      if '+' in statement:
        x += 1
      else:
        x -= 1
        
    # Print the final value of x.
    print(x)
    
  except (ValueError, IndexError):
    # Handle potential empty input or malformed numbers.
    print("Invalid input. Please provide the input as described in the problem.")

solve_bit_plus_plus()
