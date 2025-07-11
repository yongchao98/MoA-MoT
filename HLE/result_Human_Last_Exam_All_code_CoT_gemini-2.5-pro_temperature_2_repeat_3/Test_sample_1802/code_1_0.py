import sys

def calculate_max_rank(d):
  """
  Calculates the maximal rank of the complementary channel of a Pauli channel
  for a d-dimensional quantum system (qudit).

  Args:
    d: The dimension of the qudit system.

  Returns:
    The maximal rank.
  """
  if not isinstance(d, int) or d < 2:
    raise ValueError("The dimension 'd' must be an integer greater than or equal to 2.")
  
  # The maximal rank of a Pauli channel on a d-dimensional system is d^2.
  # The rank of the complementary channel is equal to the rank of the original channel.
  rank = d**2
  return rank

def main():
  """
  Main function to execute the calculation and print the result.
  """
  try:
    # Set the dimension 'd' for the qudit.
    # You can change this value to calculate the rank for a different dimension.
    d = 3
    
    # Or take 'd' from command line arguments for more interactivity
    if len(sys.argv) > 1:
        try:
            d = int(sys.argv[1])
        except ValueError:
            print(f"Invalid input '{sys.argv[1]}'. Using default d={d}.")


    max_rank = calculate_max_rank(d)
    
    print("For a qudit Pauli channel on a system of dimension d:")
    print(f"d = {d}")
    print("\nThe maximal rank of the complementary channel is d^2.")
    print(f"Calculation: {d} * {d} = {max_rank}")

  except (ValueError, TypeError) as e:
    print(f"Error: {e}", file=sys.stderr)
    
if __name__ == "__main__":
  main()