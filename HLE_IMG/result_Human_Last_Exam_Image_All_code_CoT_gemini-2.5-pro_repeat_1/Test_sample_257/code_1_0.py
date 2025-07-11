def analyze_nmr_peak():
    """
    Analyzes the 1H NMR spectrum of Compound 1 based on its predicted structure.

    The reaction of Pr-DAOTA with concentrated sulfuric acid leads to a symmetrical
    disulfonated product. This analysis focuses on the most deshielded proton in that product.
    """

    # Step 1: Identify the most deshielded proton and count its neighbors.
    # The most deshielded proton is the lone proton on the central, electron-deficient
    # pyridinium-like ring. It has no protons on adjacent carbons.
    n = 0  # Number of adjacent, non-equivalent protons

    # Step 2: Calculate the splitting pattern using the n+1 rule.
    # The multiplicity (number of peaks) is given by the equation: Multiplicity = n + 1
    multiplicity = n + 1

    if multiplicity == 1:
        splitting_pattern = "singlet (s)"
    elif multiplicity == 2:
        splitting_pattern = "doublet (d)"
    elif multiplicity == 3:
        splitting_pattern = "triplet (t)"
    elif multiplicity == 4:
        splitting_pattern = "quartet (q)"
    else:
        splitting_pattern = "multiplet (m)"

    # Step 3: Determine the integration.
    # There is only one such proton in the entire molecule.
    integration = 1

    # Step 4: Print the final answer, showing the logic.
    print("Analysis of the highest deshielded proton peak in Compound 1:")
    print("-" * 50)

    print("Splitting Pattern Calculation (based on the n+1 rule):")
    print(f"The number of adjacent protons (n) is: {n}")
    print(f"The multiplicity is calculated as n + 1 = {n} + 1 = {multiplicity}")
    print(f"A multiplicity of {multiplicity} corresponds to a: {splitting_pattern}")

    print("\nIntegration Determination:")
    print(f"The number of protons represented by this peak is: {integration}")
    print("-" * 50)

    print(f"\nFinal Conclusion: The splitting pattern is a {splitting_pattern} and the integration is {integration}H.")


if __name__ == "__main__":
    analyze_nmr_peak()