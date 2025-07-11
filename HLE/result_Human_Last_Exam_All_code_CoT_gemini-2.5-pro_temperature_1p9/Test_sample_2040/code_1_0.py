import textwrap

def analyze_spdc_potential(materials):
    """
    Analyzes the potential for materials to exhibit Spontaneous Parametric
    Down-Conversion (SPDC) based on their crystal symmetry.
    """
    print("--- Analysis of Spontaneous Parametric Down-Conversion (SPDC) in Boron Nanosheets ---\n")

    # SPDC is a second-order nonlinear optical process.
    # A fundamental requirement for a material to exhibit this effect
    # is that it must lack a center of inversion symmetry.
    # Materials with inversion symmetry are 'centrosymmetric'.
    # Centrosymmetric materials have a second-order nonlinear susceptibility (χ⁽²⁾) of zero.
    spdc_rule = "SPDC requires a non-centrosymmetric crystal structure (i.e., inversion_symmetry = False)."
    print(f"Physical Principle: {spdc_rule}\n")

    spdc_expected_overall = False

    print("Investigating common free-standing borophene phases:")
    for key, props in materials.items():
        name = props['name']
        space_group = props['space_group']
        has_inversion_symmetry = props['inversion_symmetry']

        # Determine if SPDC is expected for this phase
        if not has_inversion_symmetry:
            spdc_expected_overall = True
            conclusion = f"is NOT centrosymmetric. Therefore, it IS expected to exhibit SPDC."
        else:
            conclusion = f"IS centrosymmetric. Therefore, it is NOT expected to exhibit SPDC."

        print(f"- {name} (space group: {space_group}):")
        print(f"  - Inversion Symmetry: {has_inversion_symmetry}")
        print(f"  - Conclusion: This structure {conclusion}\n")

    print("--- Overall Conclusion ---")
    if not spdc_expected_overall:
        main_conclusion = (
            "Based on the crystal symmetry of its common stable phases, "
            "ideal, free-standing boron nanosheets would NOT be expected to exhibit "
            "spontaneous parametric down-conversion."
        )
        print(main_conclusion)
    else:
        # This case is not reached with the current data but is included for completeness
        main_conclusion = (
            "Some phases of borophene are non-centrosymmetric and would be expected to exhibit SPDC."
        )
        print(main_conclusion)

    print("\n--- Important Caveats ---")
    caveats = (
        "This analysis is based on the ideal, defect-free crystal structure. In real materials, "
        "factors such as strain, wrinkles, edges, vacancies, or the formation of other "
        "non-centrosymmetric polymorphs could break the inversion symmetry. This breakage could "
        "potentially induce a weak second-order nonlinear response, including SPDC, even in "
        "theoretically centrosymmetric materials."
    )
    # Using textwrap for clean printing in terminals
    print(textwrap.fill(caveats, width=80))

# --- Data Section ---
# Properties for common free-standing borophene phases.
# Research shows the most stable predicted structures (β12 and χ3) are both centrosymmetric.
material_database = {
    'borophene_beta12': {
        'name': 'β12 Borophene',
        'space_group': 'Pmmn',
        'inversion_symmetry': True
    },
    'borophene_chi3': {
        'name': 'χ3 Borophene',
        'space_group': 'Cmcm',
        'inversion_symmetry': True
    }
}

# Run the analysis
if __name__ == "__main__":
    analyze_spdc_potential(material_database)