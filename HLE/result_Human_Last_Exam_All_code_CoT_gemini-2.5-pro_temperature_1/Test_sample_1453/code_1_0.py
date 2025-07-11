import textwrap

def analyze_coiled_coil(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.
    """
    print("Step 1: The protein sequence to be analyzed is:")
    print(f"'{sequence}'\n")

    print("Step 2: Understanding the Heptad Repeat.")
    explanation = (
        "Coiled-coils are characterized by a repeating seven-amino-acid pattern "
        "called a heptad repeat, labeled (abcdefg). The 'a' and 'd' positions "
        "are typically hydrophobic and form the core of the structure. The "
        "specific amino acids at these positions determine the oligomeric state."
    )
    print(textwrap.fill(explanation, width=80))
    print("\nStep 3: Testing all 7 possible registers to find the correct alignment.")
    print("-" * 60)

    possible_registers = []
    for start_pos_of_a in range(1, 8):
        a_residues = []
        d_residues = []
        # In 0-based indexing, position 'a' is at index (start_pos_of_a - 1)
        # and moves in steps of 7.
        # Position 'd' is 3 positions after 'a'.
        for i in range(len(sequence)):
            # The position within the heptad (0-6 for a-g)
            heptad_pos = (i - (start_pos_of_a - 1)) % 7
            if heptad_pos == 0:  # 'a' position
                a_residues.append(sequence[i])
            elif heptad_pos == 3:  # 'd' position
                d_residues.append(sequence[i])

        print(f"Register where residue 1 is position '{'abcdefg'[-(start_pos_of_a-1)]}':")
        print(f"  'a' positions have residues: {', '.join(a_residues)}")
        print(f"  'd' positions have residues: {', '.join(d_residues)}\n")
        possible_registers.append({'a': a_residues, 'd': d_residues})
    print("-" * 60)

    print("Step 4: Interpretation of the results.")
    interpretation_1 = (
        "We are looking for a register where both 'a' and 'd' positions are "
        "consistently filled with hydrophobic residues. Looking at the results "
        "above, the register where residue 1 is 'e' (and thus residue 4, 'A', is 'a') "
        "is the most plausible."
    )
    print(textwrap.fill(interpretation_1, width=80))
    print("\nFor this optimal register:")
    print("  'a' positions are consistently Alanine (A).")
    print("  'd' positions are consistently Leucine (L).\n")

    interpretation_2 = (
        "This 'a=Alanine, d=Leucine' pattern is a classic indicator of a trimeric "
        "(3-helix) coiled-coil. The small side chain of Alanine at the 'a' "
        "position allows the helices to pack together more closely than in a dimer, "
        "creating a perfectly packed hydrophobic core with the large Leucine "
        "residues at the 'd' position. This arrangement follows the "
        "'knobs-into-holes' packing rule with maximum efficiency for a trimer."
    )
    print(textwrap.fill(interpretation_2, width=80))

    final_answer = 3
    print("\nStep 5: Conclusion.")
    print(f"Based on the analysis, the predicted oligomeric state is {final_answer}.")


# The protein sequence in question
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

analyze_coiled_coil(protein_sequence)

<<<B>>>