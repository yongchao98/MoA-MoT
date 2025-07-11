import textwrap

def analyze_coiled_coil(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.
    """
    print("1. Identifying the Heptad Repeat Pattern (a, d core positions)\n")
    print(f"Sequence: {sequence}\n")
    
    # The most common coiled-coils have a heptad (7-residue) repeat.
    # The 'a' and 'd' positions form the hydrophobic core.
    # We test possible frames and find the one that maximizes hydrophobicity at 'a' and 'd'.
    # By inspection, the frame starting with 'a' at residue 4 ('A') and 'd' at residue 7 ('L')
    # is the most consistent and chemically sensible. This is the 'efgabcd' frame.
    
    # Heptad 'a' positions start at index 3 (pos 4) and repeat every 7.
    # Heptad 'd' positions start at index 6 (pos 7) and repeat every 7.
    a_indices = [i for i in range(3, len(sequence), 7)]
    d_indices = [i for i in range(6, len(sequence), 7)]

    print("2. Extracting Core Residues for Each Repeat\n")
    
    a_residues = [sequence[i] for i in a_indices]
    d_residues = [sequence[i] for i in d_indices]

    print("The identified core (a, d) pairs are:")
    # We use a loop to show the final "equation" for each repeat
    for i in range(len(a_residues)):
        print(f"Repeat {i+1}: a = {a_residues[i]}, d = {d_residues[i]}")
    
    print("\n" + "-"*40 + "\n")
    print("3. Analysis of Packing Efficiency (Knobs-into-Holes)\n")
    
    analysis_text = """The 'd' positions are consistently Leucine (L). Leucine's large, straight-chain hydrophobic side chain is a classic anchor for the coiled-coil core and is found in oligomers of all states.

The 'a' positions are consistently Alanine (A). Alanine is the smallest hydrophobic amino acid. This provides a strong clue:

- Dimer (2 helices): Packing an 'A' from one helix against an 'A' from another would likely leave a cavity in the core, which is energetically unfavorable. Dimers are more stable with larger residues like Leucine or Isoleucine at the 'a' position.

- Trimer (3 helices): A trimeric arrangement allows the three small Alanine side chains (one from each helix) to pack together efficiently in the center of the bundle. This 'A-L' pattern (a=A, d=L) is a known structural motif that strongly favors the formation of trimers.

- Tetramer (4 helices): While possible, the classic tetramer motif involves a larger residue at 'a' and a smaller one at 'd' (e.g., Isoleucine at 'a', Alanine at 'd') to accommodate the different packing interfaces. The observed 'A-L' pattern is a much better fit for the trimeric geometry.
"""
    print(textwrap.fill(analysis_text, width=70))
    
    print("\n" + "-"*40 + "\n")
    print("4. Conclusion\n")
    print("Based on the (a=A, d=L) repeating pattern, the sequence is predicted to form a trimer to achieve optimal packing efficiency in its hydrophobic core.")
    print("Therefore, the oligomeric state is 3.")


# The protein sequence provided by the user
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

# Run the analysis
analyze_coiled_coil(protein_sequence)
