def explain_helix_type():
    """
    Explains the reasoning for identifying the helix type of an alternating alpha/epsilon-peptidomimetic.
    """
    # 1. Description of the molecule and the problem.
    print("Step 1: Analyzing the Peptidomimetic's Composition")
    print("The molecule is a foldamer, an artificial peptide, with an alternating sequence of:")
    print("  - Alanine: A standard, short 'alpha' (α) amino acid.")
    print("  - Cyclically-constrained 'epsilon' (ε) amino acid: A long-chain non-natural amino acid.")
    print("The goal is to determine the most likely helical structure it will form.\n")

    # 2. Principles of hybrid peptide folding.
    print("Step 2: Principles of Foldamer Helices")
    print("Helices are stabilized by hydrogen bonds (H-bonds) that create rings of atoms.")
    print("The type of helix is defined by the number of atoms in these rings.")
    print("Alternating copolymers like this often form helices with two different H-bond ring sizes, denoted as 'X/Y'.")
    print("The size of the rings is directly related to the length of the amino acid building blocks.\n")

    # 3. Relating residue size to ring size.
    print("Step 3: Predicting Ring Size from Residue Length")
    print("As the amino acid backbone gets longer, the H-bonded rings get larger:")
    print("  - α-peptides -> ~13-membered rings")
    print("  - β-peptides -> ~14-membered rings")
    print("  - γ-peptides -> ~16-membered rings")
    print("  - δ-peptides -> ~18-membered rings")
    print("Since our molecule contains a very long ε-residue, we expect very large rings.\n")
    
    # 4. Final deduction and explanation of the 18/20-helix
    print("Step 4: Conclusion based on Evidence")
    print("Among the choices, '18/20' represents a helix with the largest rings, consistent with the presence of ε-amino acids.")
    print("This is strongly supported by experimental studies in chemical literature, which have established that alternating α/ε-peptides form a stable '18/20-helix'.\n")

    # 5. Output the final "equation" describing the bonding pattern
    print("The final 'equation' describes the two types of hydrogen bonds forming the helix:")
    print("-" * 70)
    # The final equation as requested.
    print("Helix Type: 18/20-Helix")
    print("This is defined by two alternating hydrogen bonds:")
    print("1. C=O group of an α-residue(i) H-bonds with the N-H group of an ε-residue(i+3) = 18-membered ring")
    print("2. C=O group of an ε-residue(i) H-bonds with the N-H group of an α-residue(i+3) = 20-membered ring")
    print("-" * 70)

explain_helix_type()