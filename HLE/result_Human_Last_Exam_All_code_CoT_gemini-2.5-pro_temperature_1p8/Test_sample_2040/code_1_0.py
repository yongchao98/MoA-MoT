import textwrap

def explain_spdc_in_borophene():
    """
    Analyzes the likelihood of Spontaneous Parametric Downconversion (SPDC)
    in free-standing boron nanosheets based on their crystal symmetry.
    """

    # Helper function for printing formatted text
    def print_step(title, content):
        print(f"--- {title} ---")
        print(textwrap.fill(content, width=80))
        print()

    # Step 1: Explain the fundamental requirement for SPDC
    step1_title = "Requirement for Spontaneous Parametric Downconversion (SPDC)"
    step1_content = (
        "SPDC is a second-order nonlinear optical process. A material can only "
        "exhibit this process if it has a non-zero second-order nonlinear "
        "susceptibility, commonly denoted as χ⁽²⁾."
    )
    print_step(step1_title, step1_content)

    # Step 2: Link the requirement to crystal symmetry
    step2_title = "The Role of Crystal Symmetry"
    step2_content = (
        "For the bulk electric dipole contribution to χ⁽²⁾ to be non-zero, "
        "the material's crystal lattice must lack a center of inversion symmetry. "
        "Materials that have a center of inversion are called 'centrosymmetric', "
        "and for these materials, χ⁽²⁾ is zero. Consequently, SPDC is forbidden "
        "in centrosymmetric materials."
    )
    print_step(step2_title, step2_content)

    # Step 3: Analyze the symmetry of boron nanosheets
    step3_title = "Symmetry of Free-Standing Boron Nanosheets (Borophene)"
    step3_content = (
        "The most stable and commonly synthesized phases of free-standing boron "
        "nanosheets, such as the β₁₂ and χ₃ phases, have been shown to be "
        "centrosymmetric. Their crystal structures (space group Pmmn) possess "
        "a center of inversion symmetry."
    )
    print_step(step3_title, step3_content)

    # Step 4: Draw the final conclusion
    step4_title = "Conclusion"
    step4_content = (
        "Because the stable phases of free-standing boron nanosheets are "
        "centrosymmetric, their second-order nonlinear susceptibility (χ⁽²⁾) is zero. "
        "Therefore, they would not be expected to exhibit spontaneous parametric "
        "downconversion. (Note: This applies to the ideal, defect-free sheet. "
        "Symmetry breaking at edges, due to strain, or interaction with a "
        "substrate could potentially induce a weak second-order response.)"
    )
    print_step(step4_title, step4_content)

if __name__ == '__main__':
    explain_spdc_in_borophene()