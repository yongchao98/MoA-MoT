def solve_nmr_puzzle():
    """
    This function determines and prints the splitting pattern and integration
    of the most deshielded proton in Compound 1.

    Reasoning:
    1.  The reaction is an electrophilic sulfonation, adding an -SO3H group to an outer ring.
    2.  The most deshielded proton is the lone proton on the central, electron-deficient pyridine-like ring.
    3.  Integration: There is only one such proton in the molecule. Therefore, the integration is 1H.
    4.  Splitting Pattern: The adjacent carbons are bridgehead carbons with no protons.
        - Number of neighboring protons (n) = 0.
        - According to the n+1 rule, the splitting pattern is n+1 = 0+1 = 1, which is a singlet.
    """
    
    # Properties of the most deshielded proton peak
    splitting_pattern = "singlet"
    integration_value = 1  # Represents 1H

    print("Analysis of the highest deshielded proton peak in Compound 1:")
    print(f"Splitting Pattern: {splitting_pattern}")
    print(f"Integration: {integration_value}H")

solve_nmr_puzzle()