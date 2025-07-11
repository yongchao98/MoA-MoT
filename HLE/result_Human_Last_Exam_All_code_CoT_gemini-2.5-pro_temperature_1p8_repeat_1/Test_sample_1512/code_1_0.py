import sys

def solve():
  """
  This function reads two large numbers as strings from standard input,
  multiplies them, and prints the result.
  Python's built-in support for arbitrary-precision integers makes this
  calculation straightforward.
  """
  try:
    # Read the first number p
    p_str = sys.stdin.readline().strip()
    # Read the second number q
    q_str = sys.stdin.readline().strip()

    # Convert the string representations to integers
    p = int(p_str)
    q = int(q_str)

    # Multiply the two integers
    o = p * q
    
    # Print the result 'o' as specified in the problem's Output section
    print(o)

  except (IOError, ValueError) as e:
    # Handle potential input errors
    print(f"An error occurred: {e}", file=sys.stderr)

solve()