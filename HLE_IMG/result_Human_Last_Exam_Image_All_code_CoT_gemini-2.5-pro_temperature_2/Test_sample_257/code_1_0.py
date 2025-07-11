def analyze_nmr_signal():
  """
  This script analyzes the 1H NMR signal for the most deshielded proton in Compound 1.

  Step 1: The most deshielded protons are the two on the central pyridine ring due to the cationic charge.
  Step 2: Due to C2 symmetry, these two protons are chemically equivalent.
  Step 3: The integration value is therefore 2.
  Step 4: The splitting pattern is determined by the number of neighbors (n).
           Each of these protons has one neighbor via a 4-bond peri-coupling. So, n=1.
  Step 5: Using the n+1 rule for splitting.
  """

  integration = 2
  number_of_neighbors = 1
  
  # Calculate splitting based on the n+1 rule
  splitting_peaks = number_of_neighbors + 1
  splitting_pattern = "doublet"

  print("The peak for the most deshielded protons has the following characteristics:")
  print(f"Integration: {integration}H")
  print("\nSplitting calculation (n+1 rule):")
  print(f"n (number of neighbors) = {number_of_neighbors}")
  print(f"Number of peaks = {number_of_neighbors} + 1 = {splitting_peaks}")
  print(f"Splitting Pattern: {splitting_pattern}")

analyze_nmr_signal()