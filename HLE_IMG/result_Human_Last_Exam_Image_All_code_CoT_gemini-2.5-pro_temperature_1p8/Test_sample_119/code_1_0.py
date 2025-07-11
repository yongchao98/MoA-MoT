# Based on the analysis of all provided spectra (MS, IR, 1H NMR, 13C NMR, DEPT-135, and HSQC),
# the structure of the compound has been determined.
# This script will print the final IUPAC name of the compound.

def get_iupac_name():
  """
  Returns the IUPAC name of the compound based on spectral data analysis.
  
  Analysis Summary:
  - MS: M+ = 135 (odd, suggests Nitrogen).
  - Formula: C9H13N (DBE = 4).
  - NMR Fragments: Phenyl group (C6H5), CH, CH2, CH3.
  - Connectivity: C6H5-CH2-CH(NH2)-CH3.
  - The structure is consistent with all spectra.
  """
  iupac_name = "1-phenylpropan-2-amine"
  return iupac_name

if __name__ == "__main__":
  final_name = get_iupac_name()
  print(f"The IUPAC name of the compound is: {final_name}")