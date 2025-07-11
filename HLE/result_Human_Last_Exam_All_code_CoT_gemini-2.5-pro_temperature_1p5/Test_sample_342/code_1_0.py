import math

def calculate_Dg(g):
  """
  Calculates the degree D_g based on the formula 2**(g-1) * (2**g + 1).
  """
  # Using math.pow returns float, so we cast to int for clean output.
  # Bitwise shifts are also a clean way to do powers of 2.
  term1 = 1 << (g - 1) # 2**(g-1)
  term2 = (1 << g) + 1  # 2**g + 1
  return term1 * term2

def main():
  """
  Main function to calculate and print the first 4 terms of the sequence D_g.
  """
  print("The sequence of the first 4 terms of D_g is derived from the formula: D_g = 2^(g-1) * (2^g + 1)")
  print("The calculation for each term is as follows:")

  g_values = range(1, 5)
  results = []
  for g in g_values:
    # We calculate the components of the formula to display them.
    part1_val = 1 << (g-1)
    power_of_2_g = 1 << g
    part2_val = power_of_2_g + 1
    result = part1_val * part2_val
    results.append(result)
    
    # Print the equation with all the numbers.
    print(f"For g={g}: D_{g} = 2^({g}-1) * (2^{g} + 1) = {part1_val} * ({power_of_2_g} + 1) = {result}")

  # Join the final numbers to present the sequence.
  final_sequence = ", ".join(map(str, results))
  print(f"\nThe first 4 terms of the sequence D_g are: {final_sequence}")

if __name__ == "__main__":
  main()