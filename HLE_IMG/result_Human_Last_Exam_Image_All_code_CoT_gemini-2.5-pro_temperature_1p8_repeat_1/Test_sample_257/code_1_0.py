def solve_nmr_puzzle():
    """
    This function determines the splitting pattern and integration
    of the most deshielded proton in Compound 1.
    """
    # Based on the chemical structure and reaction, the most deshielded proton is
    # the unique proton on the central aromatic ring of the DAOTA core.
    
    # Integration is determined by the number of protons giving rise to the signal.
    # There is only one such unique proton in the molecule.
    integration = 1
    
    # The splitting pattern is determined by the number of adjacent equivalent protons (n)
    # using the n+1 rule. This proton has two equivalent neighbors.
    num_neighbors = 2
    multiplicity_value = num_neighbors + 1
    
    if multiplicity_value == 3:
        splitting_pattern = "triplet"
    else:
        # Fallback for other cases, though triplet is correct here.
        splitting_pattern = "multiplet"
        
    print(f"The splitting pattern of the highest deshielded proton is a {splitting_pattern}.")
    print(f"The integration of this peak is {integration}H.")

solve_nmr_puzzle()