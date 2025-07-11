def get_nmr_properties():
    """
    This function determines and prints the splitting pattern and integration
    of the most deshielded proton in Compound 1.

    - The most deshielded proton is the lone proton on the central, electron-deficient
      aromatic ring of the diazaoxatriangulenium core.
    - Splitting Pattern: It has no adjacent protons (n=0), so by the n+1 rule,
      the splitting pattern is a singlet (0+1=1).
    - Integration: There is only one such proton in the molecule, so the
      integration is 1.
    """
    
    # Number of neighboring protons
    n = 0
    
    # Calculate the number of peaks in the signal (n+1 rule)
    number_of_peaks = n + 1
    
    if number_of_peaks == 1:
        splitting_pattern = "singlet"
    elif number_of_peaks == 2:
        splitting_pattern = "doublet"
    elif number_of_peaks == 3:
        splitting_pattern = "triplet"
    else:
        splitting_pattern = "multiplet"
        
    # Number of protons corresponding to the signal
    integration = 1
    
    print(f"The splitting pattern of the highest deshielded proton peak is a {splitting_pattern}.")
    print(f"The integration of this peak is {integration}H.")

get_nmr_properties()