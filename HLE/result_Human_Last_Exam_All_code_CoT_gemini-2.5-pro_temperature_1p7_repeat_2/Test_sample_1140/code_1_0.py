import sys

def synthesis_plan():
    """
    This function explains and outlines the Native Chemical Ligation (NCL)
    strategy for synthesizing a large peptide with an unnatural amino acid.
    """
    total_length = 100
    # NCL requires a Cysteine residue at the ligation point.
    # We will assume a Cysteine is strategically placed around the midpoint.
    ligation_site_cys_position = 45

    # Calculate the end and start points for the two fragments
    fragment_1_end = ligation_site_cys_position - 1
    fragment_2_start = ligation_site_cys_position

    print("The most helpful technique for this task is Native Chemical Ligation (NCL).")
    print("\nREASONING:")
    print("A 100 amino acid peptide is too long for efficient direct synthesis by standard Solid-Phase Peptide Synthesis (SPPS).")
    print("NCL circumvents this by chemically joining two smaller, separately synthesized peptide fragments. This method is highly compatible with unnatural amino acids.")
    print("\nThe provided sequence snippet contains a Cysteine (C), which is the required site for the NCL reaction.")

    print("\n--- SYNTHESIS PLAN ---")
    print(f"1. Synthesize Peptide Fragment 1 (residues 1-{fragment_1_end}) via SPPS. This fragment is modified to have a C-terminal thioester.")
    print(f"2. Synthesize Peptide Fragment 2 (residues {fragment_2_start}-{total_length}) via SPPS. This fragment starts with the native Cysteine and will contain the unnatural amino acid (Azido Phenylalanine).")
    print("3. Purify both fragments to high homogeneity.")
    print("4. Ligate Fragment 1 and Fragment 2 in solution. The N-terminal Cysteine of Fragment 2 attacks the C-terminal thioester of Fragment 1, forming the final full-length peptide.")

    print("\n--- LIGATION EQUATION ---")
    # The prompt requests printing each number in the final equation.
    # We will build the string piece by piece to demonstrate this.
    sys.stdout.write("Peptide(")
    sys.stdout.write(str(1))
    sys.stdout.write("-")
    sys.stdout.write(str(fragment_1_end))
    sys.stdout.write(")-Thioester + Cys-Peptide(")
    sys.stdout.write(str(fragment_2_start))
    sys.stdout.write("-")
    sys.stdout.write(str(total_length))
    sys.stdout.write(")  --->  Full-Peptide(")
    sys.stdout.write(str(1))
    sys.stdout.write("-")
    sys.stdout.write(str(total_length))
    sys.stdout.write(")\n")


# Execute the function to print the plan
synthesis_plan()