import sys

def explain_spdc_in_boron_nanosheets():
    """
    This script explains, based on physical principles, whether free-standing boron
    nanosheets would be expected to exhibit spontaneous parametric down-conversion (SPDC).
    """

    # Introduction to the core physical process
    print("Analysis: Can Boron Nanosheets Exhibit Spontaneous Parametric Down-Conversion (SPDC)?")
    print("-" * 80)
    print("1. SPDC is a second-order nonlinear optical process.")
    print("The optical polarization (P) response of a material to an intense electric field (E) is described by the equation:")
    print("   P = ε₀ * (χ(1)E + χ(2)E² + χ(3)E³ + ...)")
    print("\nFor SPDC to occur, the material's second-order susceptibility, χ(2), must be non-zero.")
    print("-" * 80)

    # The role of crystal symmetry
    print("2. The fundamental requirement is a lack of inversion symmetry.")
    print("In materials with a center of symmetry (centrosymmetric materials), the polarization must reverse perfectly if the electric field is reversed (P(-E) = -P(E)).")
    print("The even-powered terms in the polarization equation, like χ(2)E², do not change sign when E -> -E.")
    print("To satisfy the symmetry condition P(-E) = -P(E), all even-order coefficients must be zero.")
    
    # State the value of chi-2 for centrosymmetric materials
    chi_2_value = 0
    print("\nTherefore, for any centrosymmetric material, the following must be true:")
    print(f"   χ(2) = {chi_2_value}")
    print("-" * 80)

    # Analysis of the specific material
    print("3. Examining the structure of free-standing boron nanosheets (borophene).")
    print("Borophene is a 2D material with several known stable crystal structures (polymorphs).")
    print("The most commonly studied phases, such as the β₁₂ and χ₃ phases, are planar and possess a center of inversion.")
    print("This means their crystal structures are centrosymmetric.")
    print("-" * 80)

    # Final conclusion
    print("4. Conclusion")
    print("Since the common phases of free-standing boron nanosheets are centrosymmetric, their second-order susceptibility (χ(2)) is zero by symmetry.")
    print("As a result, they are NOT expected to exhibit second-order phenomena like spontaneous parametric downconversion.")

# Execute the explanation
explain_spdc_in_boron_nanosheets()