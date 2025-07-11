def solve_nmr_puzzle():
    """
    This script determines the splitting pattern and integration for the most
    deshielded proton in Compound 1 based on chemical principles.
    
    The reaction is the sulfonation of Pr-DAOTA. The reaction will occur
    on the electron-rich outer rings, while the electron-deficient central
    ring remains intact. The most deshielded proton will be on this central
    positively charged ring. Due to the molecule's symmetry, this unique
    proton (let's call it Hd) is flanked by two other equivalent protons (Hc).
    """

    # Identify the most deshielded proton and its neighbors
    # The most deshielded proton is the unique one on the central aromatic ring.
    # It has 2 equivalent neighboring protons.
    num_adjacent_protons = 2

    # Calculate the splitting pattern using the n+1 rule
    # multiplicity = n + 1
    splitting_multiplicity = num_adjacent_protons + 1
    
    if splitting_multiplicity == 1:
        splitting_pattern_name = "singlet"
    elif splitting_multiplicity == 2:
        splitting_pattern_name = "doublet"
    elif splitting_multiplicity == 3:
        splitting_pattern_name = "triplet"
    elif splitting_multiplicity == 4:
        splitting_pattern_name = "quartet"
    else:
        splitting_pattern_name = "multiplet"

    # Determine the integration
    # There is only one such proton in the entire molecule.
    integration_value = 1

    # Print the step-by-step conclusion
    print("Step-by-step determination:")
    print(f"1. The most deshielded proton is the one on the central acridinium-like ring.")
    print(f"2. This proton has {num_adjacent_protons} equivalent neighboring protons.")
    print(f"3. According to the n+1 rule, its splitting pattern is (n+1) = {num_adjacent_protons} + 1 = {splitting_multiplicity}, which is a {splitting_pattern_name}.")
    print(f"4. As there is only one such proton in the molecule, its integration is {integration_value}H.")
    
    # Print the final answer clearly
    print("\n--- Final Answer ---")
    print(f"The splitting pattern is a {splitting_pattern_name}.")
    print(f"The integration is {integration_value}H.")


solve_nmr_puzzle()