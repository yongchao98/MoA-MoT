def analyze_nmr_peak():
    """
    This script analyzes the provided reaction to determine the splitting
    pattern and integration of the most deshielded proton in the product, Compound 1.
    """

    print("--- Analysis of the Most Deshielded Proton in Compound 1 ---")

    # Step 1: Identify the most deshielded proton in the molecule.
    # The most deshielded proton is the one in the most electron-poor environment.
    # In Compound 1, this is the single proton on the central aromatic ring,
    # which is part of a positively charged, electron-deficient system.
    most_deshielded_proton_identity = "The single proton on the central aromatic ring"
    print(f"1. Identity of the most deshielded proton: {most_deshielded_proton_identity}")

    # Step 2: Determine the splitting pattern using the n+1 rule.
    # Splitting is caused by coupling to protons on adjacent atoms.
    # This proton has no protons on its neighboring carbons.
    neighboring_protons_n = 0
    print(f"\n2. Determining the Splitting Pattern:")
    print(f"   - Number of neighboring protons (n) = {neighboring_protons_n}")

    # Calculate the multiplicity based on the n+1 rule.
    multiplicity = neighboring_protons_n + 1
    print(f"   - Applying the n+1 rule: Multiplicity = n + 1 = {neighboring_protons_n} + 1 = {multiplicity}")

    if multiplicity == 1:
        splitting_pattern = "singlet"
    elif multiplicity == 2:
        splitting_pattern = "doublet"
    elif multiplicity == 3:
        splitting_pattern = "triplet"
    else:
        splitting_pattern = "multiplet"

    print(f"   - The calculated splitting pattern is a '{splitting_pattern}'.")

    # Step 3: Determine the integration of the peak.
    # Integration represents the number of protons giving rise to the signal.
    # By chemical structure, there is only one proton of this type.
    integration_value = 1
    print(f"\n3. Determining the Integration:")
    print(f"   - Number of protons of this type in the molecule = {integration_value}")
    print(f"   - The integration for this peak is {integration_value}H.")

    # Step 4: Final conclusion.
    print("\n--- Conclusion ---")
    print("The highest deshielded proton peak in the 1H NMR spectrum of Compound 1 is a "
          f"{splitting_pattern} with an integration of {integration_value}H.")


if __name__ == '__main__':
    analyze_nmr_peak()