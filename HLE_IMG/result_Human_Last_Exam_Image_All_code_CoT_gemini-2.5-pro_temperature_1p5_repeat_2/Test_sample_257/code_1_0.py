def solve_nmr_puzzle():
    """
    This script determines the splitting pattern and integration for the most
    deshielded proton in Compound 1 based on its molecular structure.

    1.  The most deshielded proton is the one on the central aromatic ring,
        between the two nitrogen atoms, due to its electron-deficient environment.
    2.  Integration is the count of protons for that signal. There is only one such proton.
    3.  Splitting is determined by the n+1 rule, where n is the number of adjacent protons.
        This proton has no adjacent protons.
    """

    # Integration corresponds to the number of protons generating the signal.
    integration = 1

    # Number of adjacent protons (n) for the most deshielded proton.
    # The adjacent carbons are bonded to nitrogen atoms and have no protons.
    n_adjacent_protons = 0

    # Calculate the multiplicity (number of peaks in the signal) using the n+1 rule.
    multiplicity = n_adjacent_protons + 1

    if multiplicity == 1:
        splitting_pattern = "singlet"
    elif multiplicity == 2:
        splitting_pattern = "doublet"
    elif multiplicity == 3:
        splitting_pattern = "triplet"
    else:
        splitting_pattern = "multiplet"

    print("Analysis of the highest deshielded proton peak in Compound 1:")
    print(f"Integration = {integration}H")
    print(f"Splitting Pattern = {splitting_pattern} (since n={n_adjacent_protons}, multiplicity = {n_adjacent_protons} + 1 = {multiplicity})")
    print("\nTherefore, the splitting pattern and integration are a singlet integrating to 1H.")

solve_nmr_puzzle()