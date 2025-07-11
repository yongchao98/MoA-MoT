import textwrap

def print_design_principles():
    """
    This script outlines the recommended molecular features for designing a
    single-ring liquid crystal with room temperature transitions.
    """
    title = "--- Plan for Designing a Room-Temperature Single-Ring Liquid Crystal ---"
    print(title)
    print("-" * len(title))

    # Summarize the goal based on the user's request
    print("\n[1. Design Objective]")
    objective_text = """
    To design a molecule that has a single benzene ring, exhibits liquid crystal
    phases (nematic or smectic), and has phase transition temperatures around
    room temperature (approximately 20 to 25 degrees Celsius).
    """
    print(textwrap.fill(objective_text, width=70))

    # Detail the recommended molecular features, which is the core of the design strategy
    print("\n[2. Recommended Molecular Features]")
    features_intro = """
    To create a rod-like molecule with the desired properties, specific
    functional groups should be attached to a central benzene ring at the
    para-positions (positions 1 and 4).
    """
    print(textwrap.fill(features_intro, width=70))

    print("\n  a) Add a Flexible Alkyl Chain:")
    alkyl_text = """
    Attach a flexible alkyl chain (general formula C(n)H(2n+1)) to one end.
    This chain promotes fluidity and is crucial for the formation of the
    liquid crystal phase. For the formula C(n)H(2n+1), the numbers are n, 2, and 1.
    A good starting point for 'n' is typically between 4 and 8.
    """
    print(textwrap.fill(alkyl_text, width=70, initial_indent='    ', subsequent_indent='    '))

    print("\n  b) Add a Polar Group:")
    polar_text = """
    Attach a polar group to the other end to increase the molecule's
    anisotropy, which helps the molecules align. Common choices include:
    - Cyano (-CN)
    - Fluoro (-F)
    - Methoxy (-OCH3)
    """
    print(textwrap.fill(polar_text, width=70, initial_indent='    ', subsequent_indent='    '))

    print("\n  c) Consider Lateral Substituents for Fine-Tuning:")
    lateral_text = """
    To precisely adjust the transition temperatures, add a small substituent
    (like -F or -CH3) to the side of the benzene ring. These groups can disrupt
    packing, often lowering the melting point into the desired range.
    """
    print(textwrap.fill(lateral_text, width=70, initial_indent='    ', subsequent_indent='    '))

    print("\n" + "-" * len(title))
    print("Executing this plan provides a robust strategy for your material design.")

# Run the function to print the design plan
if __name__ == '__main__':
    print_design_principles()