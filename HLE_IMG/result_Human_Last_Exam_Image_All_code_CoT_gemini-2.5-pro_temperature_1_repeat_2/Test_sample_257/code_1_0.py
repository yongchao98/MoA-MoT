def solve_nmr_puzzle():
    """
    This script programmatically determines the 1H NMR splitting pattern and
    integration for the most deshielded proton in the reaction product.
    """

    print("Step 1: Analyzing the reaction and product structure.")
    print("The reaction is an electrophilic aromatic sulfonation, adding two -SO3H groups.")
    print("The most deshielded proton is the unique proton on the central aromatic ring.\n")

    # Step 2: Determine the integration of the most deshielded proton.
    # There is only one such proton in the symmetrical molecule.
    integration = 1
    print(f"Step 2: Determining the integration.")
    print(f"Number of most deshielded protons = {integration}")
    print(f"Resulting integration = {integration}H\n")

    # Step 3: Determine the splitting pattern using the n+1 rule.
    # This proton has two adjacent, chemically equivalent neighbors due to symmetry.
    n_equivalent_neighbors = 2
    multiplicity = n_equivalent_neighbors + 1
    
    if multiplicity == 3:
        pattern_name = "Triplet"
    elif multiplicity == 1:
        pattern_name = "Singlet"
    elif multiplicity == 2:
        pattern_name = "Doublet"
    elif multiplicity == 4:
        pattern_name = "Quartet"
    else:
        pattern_name = "Multiplet"

    print("Step 3: Determining the splitting pattern.")
    print(f"Number of equivalent neighboring protons (n) = {n_equivalent_neighbors}")
    print(f"Applying the n+1 rule: n + 1 = {n_equivalent_neighbors} + 1 = {multiplicity}")
    print(f"A signal with a multiplicity of {multiplicity} is a {pattern_name}.\n")
    
    print("---")
    print("Final Answer:")
    print(f"The highest deshielded proton peak has a splitting pattern of a {pattern_name} and an integration of {integration}H.")

# Execute the analysis
solve_nmr_puzzle()