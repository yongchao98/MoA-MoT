def analyze_spdc_in_borophene():
    """
    Analyzes and explains whether free-standing boron nanosheets
    are expected to exhibit Spontaneous Parametric Down-Conversion (SPDC).
    """

    # Step 1: Define the physical principle behind SPDC.
    principle = "Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process. The efficiency of this process is governed by the material's second-order nonlinear susceptibility, χ⁽²⁾."

    # Step 2: State the symmetry requirement for a non-zero χ⁽²⁾.
    symmetry_rule = "A fundamental requirement for a material to possess a non-zero χ⁽²⁾ is that its crystal structure must lack a center of inversion. Materials with a center of inversion (centrosymmetric materials) have χ⁽²⁾ = 0 by symmetry."

    # Step 3: Analyze the crystal structure of free-standing boron nanosheets.
    borophene_analysis = "The most commonly studied and stable phases of free-standing boron nanosheets (borophene), such as the β12 and χ3 phases, are known to be centrosymmetric. Their crystal lattices possess a center of inversion."

    # Step 4: Formulate the conclusion based on the analysis.
    conclusion = "Since the intrinsic structure of these common boron nanosheets is centrosymmetric, their bulk second-order nonlinear susceptibility (χ⁽²⁾) is zero. Therefore, free-standing boron nanosheets are not expected to exhibit spontaneous parametric downconversion."

    # Print the step-by-step reasoning
    print("Step-by-step analysis:")
    print("--------------------------------------------------------------------------")
    print(f"1. Physical Principle: {principle}")
    print("--------------------------------------------------------------------------")
    print(f"2. Symmetry Requirement: {symmetry_rule}")
    print("--------------------------------------------------------------------------")
    print(f"3. Material Analysis: {borophene_analysis}")
    print("--------------------------------------------------------------------------")
    print(f"4. Conclusion: {conclusion}")
    print("--------------------------------------------------------------------------")

if __name__ == '__main__':
    analyze_spdc_in_borophene()
