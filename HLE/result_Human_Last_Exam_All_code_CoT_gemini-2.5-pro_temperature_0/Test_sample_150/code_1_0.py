import sys

def solve_bit_plus_plus():
  """
  Solves the Bit++ problem by reading statements and calculating the final value of x.
  The initial value of x is 0.
  """
  try:
    # Read the number of statements from standard input.
    n = int(sys.stdin.readline())
    
    # Initialize the variable x to 0.
    x = 0
    
    # Process each statement.
    for _ in range(n):
      statement = sys.stdin.readline().strip()
      
      # The most efficient check is to look at the middle character.
      # If it's '+', it's an increment operation.
      # If it's '-', it's a decrement operation.
      if statement[1] == '+':
        x += 1
      else:
        x -= 1
        
    # Print the final value of x.
    print(x)
    
  except (IOError, ValueError) as e:
    # Handle potential input errors gracefully.
    print(f"An error occurred: {e}", file=sys.stderr)

solve_bit_plus_plus()