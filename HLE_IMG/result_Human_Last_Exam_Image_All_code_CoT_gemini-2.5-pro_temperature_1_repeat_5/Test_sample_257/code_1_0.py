def analyze_nmr_of_compound_1():
    """
    This script determines the splitting pattern and integration for the most
    deshielded proton in Compound 1, the product of the reaction of Pr-DAOTA
    with concentrated sulfuric acid.

    The analysis concludes:
    1.  The reaction is a sulfonation, adding -SO3H groups to the molecule.
    2.  The most deshielded proton is the unique hydrogen on the central ring,
        located in a sterically hindered 'bay region'.
    3.  This proton has no adjacent protons to couple with.
    4.  There is only one such proton in the molecule.
    """

    # Properties of the most deshielded proton peak
    splitting_pattern = "Singlet"
    integration = 1

    print("Analysis of the highest deshielded proton peak in Compound 1:")
    print(f"Splitting Pattern: {splitting_pattern}")
    print(f"Integration: {integration}H")

analyze_nmr_of_compound_1()