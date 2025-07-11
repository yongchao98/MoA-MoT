def solve_nmr_puzzle():
  """
  This function solves the NMR puzzle by determining the splitting pattern and
  integration of the most deshielded proton in Compound 1.

  The steps are as follows:
  1. The reaction is an electrophilic aromatic sulfonation.
  2. The most deshielded proton is the one on the central ring, between the two nitrogen atoms, due to the extreme electron-withdrawing environment.
  3. Integration: There is only one such proton in the molecule, so the integration is 1H.
  4. Splitting Pattern: This proton has no vicinal (3-bond) proton neighbors (n=0). According to the n+1 rule, its signal is a singlet (0+1=1).
  """

  splitting_pattern = "Singlet"
  integration = "1H"

  print(f"The most deshielded proton in Compound 1 is the single proton on the central aromatic ring.")
  print(f"Splitting Pattern: {splitting_pattern}")
  print(f"Integration: {integration}")

solve_nmr_puzzle()