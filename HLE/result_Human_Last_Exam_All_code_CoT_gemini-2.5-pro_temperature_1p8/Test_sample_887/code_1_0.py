import textwrap

def analyze_reaction_dynamics():
    """
    Analyzes the effect of vibrational excitation on the F + CHD3 reaction.
    """

    # --- 1. Define the Reactants and Pathways ---
    print("Chemical System Analysis:")
    print("-------------------------")
    print("Reactants: Atomic Fluorine (F) and Monodeuterated Methane (CHD3)")
    print("Action: An infrared laser specifically excites the C-H bond's vibrational stretch.")
    print("\nThere are two competing reaction pathways:")
    print("  1. H-atom abstraction: F + C-H(D3) -> HF + .CD3")
    print("  2. D-atom abstraction: F + C-D(HD2) -> DF + .CHD2\n")

    # --- 2. Explain the Key Principle: Mode-Specific Chemistry ---
    print("Core Concept: Mode-Specific Chemistry")
    print("------------------------------------")
    explanation = (
        "This experiment is a famous example of mode-specific chemistry. "
        "The core idea is that energy placed into a specific motion (a 'mode') "
        "of a molecule can selectively promote a reaction related to that motion. "
        "Here, energy is put directly into the C-H bond vibration, which is the "
        "exact motion needed to break that bond."
    )
    print("\n".join(textwrap.wrap(explanation, width=70)))
    print()

    # --- 3. Describe the Outcome ---
    print("Predicted Effects of C-H Excitation:")
    print("------------------------------------")
    print("1. Acceleration:")
    accel_explanation = (
        "Putting vibrational energy into the C-H bond effectively lowers the "
        "energy barrier for breaking it. This makes the H-atom abstraction pathway "
        "much faster, thereby accelerating the overall reaction rate."
    )
    print("   " + "\n   ".join(textwrap.wrap(accel_explanation, width=65)))

    print("\n2. Enhanced Selectivity:")
    select_explanation = (
        "Without excitation, the weaker C-H bond is already more likely to react than the "
        "stronger C-D bond (this is the kinetic isotope effect). By selectively 'energizing' "
        "only the C-H bond, this natural preference is dramatically increased. The reaction is "
        "steered strongly toward producing HF over DF."
    )
    print("   " + "\n   ".join(textwrap.wrap(select_explanation, width=65)))
    print()

    # --- 4. Final Conclusion based on Answer Choices ---
    print("Conclusion:")
    print("-----------")
    conclusion_text = (
        "The excitation both accelerates the reaction and enhances the likelihood of "
        "H-atom removal over D-atom removal. Therefore, the most complete and accurate "
        "description is that the laser excitation accelerates the reaction by enhancing the "
        "selectivity for C-H bond cleavage."
    )
    print("\n".join(textwrap.wrap(conclusion_text, width=70)))

if __name__ == '__main__':
    analyze_reaction_dynamics()