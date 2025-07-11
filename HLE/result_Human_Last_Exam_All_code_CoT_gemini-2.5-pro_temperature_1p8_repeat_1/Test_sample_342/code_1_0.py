def calculate_Dg(g):
  """
  Calculates the smallest degree D_g of a finite etale cover of A_g
  such that the polarisation is represented by a symmetric line bundle
  on the base-change of U_g to the cover.

  The formula is D_g = 2^(g-1) * (2^g + 1).
  """
  term_1 = 2**(g - 1)
  term_2 = 2**g + 1
  result = term_1 * term_2
  return term_1, term_2, result

def main():
  """
  Main function to calculate and print the first 4 terms of the sequence D_g.
  """
  print("The sequence D_g is the smallest degree of a finite etale cover of A_g for the universal polarisation to be symmetric.")
  print("The formula for D_g is 2^(g-1) * (2^g + 1).\n")
  
  sequence_terms = []
  for g in range(1, 5):
    t1, t2, dg = calculate_Dg(g)
    print(f"For g = {g}, the calculation is: D_{g} = 2^({g}-1) * (2^{g} + 1) = {t1} * {t2} = {dg}")
    sequence_terms.append(str(dg))
  
  print("\nThe sequence of the first 4 terms is: " + ", ".join(sequence_terms))

if __name__ == "__main__":
  main()