import textwrap

def explain_reaction():
    """
    Explains the effect of the specified reaction on the cerium atoms in Ce2@C80.
    """
    # Introduction to the components
    reactant1 = "Ce2@C80"
    reactant2 = "1,1,2,2-tetrakis(2,4,6-trimethylphenyl)-1,2-disilirane"
    product_type = "exohedrally functionalized fullerene"

    # Step-by-step explanation
    explanation = [
        f"1. The reaction is between the endohedral fullerene {reactant1} and the reagent {reactant2}.",
        "2. This reaction is an 'exohedral functionalization', meaning the disilirane reagent attaches to the OUTER surface of the C80 carbon cage. It does not enter the cage.",
        "3. In the original, unmodified Ce2@C80, the two cerium atoms are known to move around freely inside the highly symmetric cage.",
        "4. The addition of the bulky chemical group to the outside breaks the symmetry of the fullerene cage. This changes the distribution of electron density on the cage's inner surface.",
        "5. This electronic change creates a new potential energy landscape inside the cage, with specific low-energy locations (potential minima).",
        "6. The positively charged cerium ions are attracted to these stable, electron-rich locations.",
        "7. Research and X-ray crystallographic studies have shown that the two cerium atoms cease their random motion and become fixed at positions along the axis defined by the point of external addition.",
        "8. These two stable positions at opposite ends of this axis are referred to as the 'poles' of the functionalized fullerene molecule."
    ]

    print("Step-by-step analysis of the chemical reaction:")
    for step in explanation:
        # The textwrap module helps format the text for better readability.
        print(textwrap.fill(step, width=80))

    print("\nConclusion:")
    print("The cerium atoms are now positioned at the poles of the fullerene.")


if __name__ == "__main__":
    explain_reaction()
