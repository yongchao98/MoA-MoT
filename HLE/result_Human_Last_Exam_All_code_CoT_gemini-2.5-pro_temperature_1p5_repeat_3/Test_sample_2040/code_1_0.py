def analyze_spdc_in_boron_nanosheets():
    """
    Analyzes and explains whether free-standing boron nanosheets would
    exhibit spontaneous parametric downconversion (SPDC).
    """

    print("Step 1: Define the primary requirement for SPDC.")
    print("Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process.")
    print("A fundamental requirement for any material to exhibit a second-order nonlinear effect is that its crystal structure must lack inversion symmetry (i.e., it must be non-centrosymmetric).")
    print("-" * 75)

    print("Step 2: Examine the crystal symmetry of boron nanosheets.")
    print("Free-standing boron nanosheets, also known as borophene, can have multiple possible structures (polymorphs).")
    print("However, the most commonly synthesized and theoretically stable phases (such as β12 and χ3) possess a center of inversion, meaning they are centrosymmetric.")
    print("-" * 75)

    print("Step 3: Combine the requirement and the material property to reach a conclusion.")
    print("In materials with a centrosymmetric crystal structure, the second-order nonlinear susceptibility (χ⁽²⁾) is zero by definition.")
    print("Since SPDC is a process that depends directly on a non-zero χ⁽²⁾, ideal, free-standing boron nanosheets in their common phases are not expected to exhibit SPDC.")
    print("-" * 75)
    
    print("Important Caveat:")
    print("While an ideal, infinite sheet would not exhibit SPDC, symmetry can be broken at the edges of a nanosheet, or by defects, strain, or interaction with a substrate. These effects could potentially lead to a very weak second-order nonlinear response, but SPDC would not be considered a characteristic feature of the material itself.")


if __name__ == "__main__":
    analyze_spdc_in_boron_nanosheets()