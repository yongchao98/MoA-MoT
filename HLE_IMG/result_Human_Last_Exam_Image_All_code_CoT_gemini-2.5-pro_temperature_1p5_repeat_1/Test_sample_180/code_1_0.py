import textwrap

def identify_protein_domain():
    """
    Analyzes the conservation pattern to identify the protein domain.
    The function prints a step-by-step reasoning based on the visual pattern.
    """
    print("Step 1: Analyze the visual pattern in the image.")
    print("The image displays a bar chart representing amino acid conservation. The key feature is the repeating, periodic pattern of red bars (highly conserved) and gray bars (variable).")
    print("Specifically, we can observe that conserved residues are consistently separated by a small number of variable residues.")
    print("-" * 70)

    print("Step 2: Connect the pattern to a known structural motif.")
    print("This type of periodic conservation pattern is the hallmark of a 'heptad repeat'.")
    print("A heptad repeat is a sequence of 7 amino acids, which can be labeled (a-b-c-d-e-f-g). This repeat is the defining feature of coiled-coil structures.")
    print("\nThe 'Leucine Zipper' is a classic example of a coiled-coil domain.")
    print("In this structure, two or more alpha-helices twist around each other. The stability of this structure comes from a hydrophobic core formed by residues at specific positions in the repeat.")
    print("\nIn the heptad repeat (a-b-c-d-e-f-g):")
    print(" - Positions 'a' and 'd' are typically hydrophobic residues (like Leucine).")
    print(" - These residues face inward, forming the hydrophobic core that 'zips' the helices together.")
    print(" - Because their role is structural, the amino acids at positions 'a' and 'd' are highly conserved.")
    print("-" * 70)

    print("Step 3: Visualize the conservation pattern.")
    print("If we represent a Conserved position with 'C' and a Variable position with 'v', the pattern for one heptad repeat is:")
    heptad_pattern = "C v v C v v v"
    positions_label = "a b c d e f g"
    print(f"Position:      ({positions_label})")
    print(f"Conservation:  ({heptad_pattern})")

    print("\nA sequence of these repeats would look like:")
    sequence_of_repeats = ' '.join([heptad_pattern] * 4)
    print(sequence_of_repeats)
    print("\nThis textual pattern, showing conserved positions at the 1st and 4th spot of every 7 residues, perfectly matches the periodic rhythm of the red bars in the image.")
    print("The other options, like Zinc Fingers (conserved Cys/His) or EGF-like domains (conserved Cys for disulfide bonds), have different and distinct conservation patterns.")
    print("-" * 70)

    print("Conclusion: The image shows the conservation pattern of a Leucine zipper motif.")


identify_protein_domain()