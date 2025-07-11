def explain_spdc_in_borophene():
    """
    This function explains whether free-standing boron nanosheets (borophene)
    are expected to exhibit Spontaneous Parametric Down-Conversion (SPDC).
    """

    # Step 1: Define SPDC and its material requirement.
    print("--- Step 1: Understanding the Physics of SPDC ---")
    print("Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process.")
    print("A critical requirement for a material to exhibit second-order effects is that its crystal structure must lack a center of inversion (it must be non-centrosymmetric).")
    print("Materials with inversion symmetry (centrosymmetric) have a second-order nonlinear susceptibility (χ⁽²⁾) of zero, which prohibits SPDC.")
    print("-" * 60)

    # Step 2: Analyze the symmetry of boron nanosheets.
    print("--- Step 2: Analyzing the Structure of Boron Nanosheets ---")
    print("Boron nanosheets, or 'borophene', can exist in several different structural forms (phases).")
    print("The question therefore is: are the common, stable phases of free-standing borophene centrosymmetric?")
    print("Scientific studies have shown that the most stable and commonly synthesized phases, such as the β₁₂ and χ₃ phases, are indeed centrosymmetric.")
    print("-" * 60)

    # Step 3: Draw the final conclusion.
    print("--- Step 3: Conclusion ---")
    print("Because the common and most stable forms of free-standing boron nanosheets possess a center of inversion symmetry, their intrinsic second-order nonlinear susceptibility (χ⁽²⁾) is zero.")
    print("Therefore, they are not expected to exhibit spontaneous parametric down-conversion.")
    print("\n(Note: SPDC could theoretically occur if a specific, less common non-centrosymmetric phase were stabilized, or if the symmetry were broken by external factors like strain or defects, but it is not an expected intrinsic property.)")
    print("-" * 60)

if __name__ == '__main__':
    explain_spdc_in_borophene()