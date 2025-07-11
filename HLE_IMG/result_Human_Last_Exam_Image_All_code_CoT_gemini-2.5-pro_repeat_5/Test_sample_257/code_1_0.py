def find_nmr_features():
    """
    This function analyzes the structure of Compound 1 to determine the
    integration and splitting pattern for its most deshielded proton.
    """
    
    # Step 1: Identify the most deshielded proton in Compound 1.
    # This is the single proton on the central aromatic ring, often labeled H5.
    
    # Step 2: Determine the integration of this proton's signal.
    # There is only one such proton in the molecule.
    integration = 1
    
    # Step 3: Determine the splitting pattern using the n+1 rule.
    # The H5 proton is coupled to the two equivalent protons on the adjacent
    # rings (H1 and H8).
    n_equivalent_neighbors = 2
    
    # Calculate the multiplicity (number of peaks).
    multiplicity = n_equivalent_neighbors + 1
    
    # Determine the name of the splitting pattern.
    if multiplicity == 3:
        splitting_pattern = "triplet"
    else:
        splitting_pattern = "multiplet"

    print("Analysis for the highest deshielded proton:")
    print(f"Integration = {integration}H")
    print(f"Number of equivalent neighboring protons (n) = {n_equivalent_neighbors}")
    print(f"Applying the n+1 rule for splitting: {n_equivalent_neighbors} + 1 = {multiplicity}")
    print(f"Splitting Pattern = {splitting_pattern}")

find_nmr_features()