def check_borophene_for_spdc():
    """
    Analyzes the potential of common boron nanosheet phases for exhibiting
    Spontaneous Parametric Down-Conversion (SPDC) based on their crystal symmetry.
    """
    print("Analyzing Boron Nanosheets for Spontaneous Parametric Down-Conversion (SPDC)\n")
    print("Principle: SPDC is a second-order nonlinear optical process whose efficiency")
    print("is primarily determined by the material's second-order nonlinear susceptibility (χ⁽²⁾).")
    print("A material can only have a non-zero χ⁽²⁾ if its crystal structure lacks a center")
    print("of inversion (i.e., it is non-centrosymmetric).\n")
    print("-" * 75)

    # Data on common borophene phases, sourced from materials science literature.
    # The key property is 'is_centrosymmetric'.
    borophene_phases = {
        "β₁₂ sheet": {
            "space_group": "Pmmn",
            "is_centrosymmetric": True,
            "comment": "A stable, buckled borophene structure."
        },
        "χ₃ sheet": {
            "space_group": "Pmmn",
            "is_centrosymmetric": True,
            "comment": "Another stable phase, also observed experimentally on Ag(111)."
        },
        "α sheet (planar honeycomb)": {
            "space_group": "P6/mmm",
            "is_centrosymmetric": True,
            "comment": "A theoretical high-symmetry planar structure."
        }
    }

    # Iterate through the phases and check their symmetry
    can_exhibit_spdc = False
    for phase, properties in borophene_phases.items():
        print(f"Checking phase: {phase}")
        print(f"  - Space Group: {properties['space_group']}")
        if properties['is_centrosymmetric']:
            print("  - Symmetry: Centrosymmetric (possesses a center of inversion)")
            print("  - Implication: Bulk χ⁽²⁾ is zero. SPDC is forbidden by symmetry.")
        else:
            print("  - Symmetry: Non-centrosymmetric")
            print("  - Implication: Bulk χ⁽²⁾ can be non-zero. SPDC is possible.")
            can_exhibit_spdc = True
        print("")

    print("-" * 75)
    print("Conclusion:")
    if not can_exhibit_spdc:
        print("Based on the crystal symmetry of their common, stable phases, free-standing")
        print("boron nanosheets are centrosymmetric.")
        print("Therefore, they are not expected to exhibit spontaneous parametric downconversion.")
    else:
        # This case is not expected based on current knowledge of stable, pure borophene
        print("Some phases of borophene are non-centrosymmetric and could potentially exhibit SPDC.")
        print("However, the most commonly studied phases are centrosymmetric.")

    print("\nNote: This analysis pertains to ideal, pure boron nanosheets. Effects from extreme")
    print("strain, defects, or chemical functionalization (which could break inversion")
    print("symmetry) are not considered here.")


# Execute the analysis
if __name__ == '__main__':
    check_borophene_for_spdc()