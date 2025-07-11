def explain_spdc_in_borophene():
    """
    This function explains whether free-standing boron nanosheets would be expected
    to exhibit spontaneous parametric downconversion (SPDC).
    """
    title = "Analysis of Spontaneous Parametric Down-Conversion (SPDC) in Boron Nanosheets"
    line = "=" * len(title)

    explanation = [
        "Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process.",
        "It involves a high-energy pump photon spontaneously splitting into a pair of lower-energy photons (signal and idler) within a nonlinear material.",
        "\nA fundamental requirement for a material to facilitate a second-order process like SPDC is that it must possess a non-zero second-order nonlinear susceptibility, denoted as chi(2).",
        "\nA material's chi(2) is non-zero only if its crystal lattice structure lacks a center of inversion symmetry (i.e., it must be non-centrosymmetric). Materials with inversion symmetry (centrosymmetric) have a chi(2) of zero in their bulk form, forbidding SPDC.",
        "\nFree-standing boron nanosheets, also known as borophene, have several known structural phases. The most stable and commonly studied phases (such as the beta_12 and chi_3 phases) are known to be centrosymmetric. They possess a center of inversion symmetry.",
        "\nTherefore, because the common, stable structures of free-standing boron nanosheets are centrosymmetric, their bulk second-order susceptibility chi(2) is zero.",
        "\nConclusion: Free-standing boron nanosheets would NOT be expected to exhibit spontaneous parametric downconversion."
    ]

    print(title)
    print(line)
    for point in explanation:
        print(f"- {point}")

# Execute the function to print the explanation
explain_spdc_in_borophene()
