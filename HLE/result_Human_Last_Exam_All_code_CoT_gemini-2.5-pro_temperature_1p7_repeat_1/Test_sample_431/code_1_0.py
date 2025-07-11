def count_cysteines_in_linker():
    """
    This script calculates the number of Cysteine residues in the
    TM3-TM4 linker of the human GABAA rho 1 receptor (UniProt: P24046).
    """

    # Full protein sequence for human GABAA rho 1 receptor (UniProt: P24046)
    protein_sequence = "MKSKAGQKSAQPKAPSVLVKRVPAVSVSSQKVERPARTVTSARTIQSYNTKNVLPVSVGSPEINAALDFDLLHFFPDTNVLYTTRLTDRLLWRDGVLLYQMFSIVTFVLNLLNVSAWTNPQRVAYATAMNIWIALACSLPGFLAIAVEYAVYEVVREVEGQGPAIADDTKATGYVDIGVYTVGQFMPTHLQSMGPLQLAQELTLEEGVKYVAINIFVFAVFLALCSEYMAAEAFEWAKQACRDLLRPVSGGPPRALALRVKIALMLCLFIFSAILAFHEYFSSKKGAAGEVKKVEAEGSGDSPRVEFFRRQGRYDFCARQILCSLEVLRCSSVSETVDRHSARIIFTFAPAFVFYVATYLLCVCMLMVSSVSAAVEVPALESMALSEKKAEAKKRKETREEILKRKGRENAL"

    # From UniProt P24046, the transmembrane domains are:
    # TM3: residues 313-333
    # TM4: residues 433-453
    # Therefore, the intracellular TM3-TM4 linker spans from residue 334 to 432.

    # In Python, strings are 0-indexed, so we subtract 1 from the start position.
    # The end position in slicing is exclusive, so 432 is correct.
    start_index = 334 - 1
    end_index = 432

    # Extract the linker sequence
    tm3_tm4_linker = protein_sequence[start_index:end_index]

    # Count the number of Cysteine ('C') residues in the linker
    cysteine_count = tm3_tm4_linker.count('C')

    print(f"Protein: Human GABAA rho 1 receptor (P24046)")
    print(f"Domain: Intracellular Linker between TM3 and TM4")
    print(f"Residue range: {start_index + 1}-{end_index}")
    print(f"Linker Sequence: {tm3_tm4_linker}")
    print(f"The number of Cysteine (C) residues found is: {cysteine_count}")

# Execute the function
count_cysteines_in_linker()
