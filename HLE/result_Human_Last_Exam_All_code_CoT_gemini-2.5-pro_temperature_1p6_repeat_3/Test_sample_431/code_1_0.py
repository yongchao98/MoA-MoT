import textwrap

def count_cysteines_in_linker():
    """
    This function counts the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAA rho-1 receptor (UniProt: P24046).
    """
    # Full amino acid sequence for human GABAA rho-1 receptor (GABRR1)
    # UniProt ID: P24046
    full_sequence = (
        "MGFVRTLLLVGLLLSHSAHASTNPDPKNTYSVDIELPQEGSVVLRTNTGAFLDINSLDV"
        "RPDTFIGYCPKDVPTCTVERSFGYAYPRLSIDVSVSSVVPNKDGKRSVAAKNWKLPKVVY"
        "ANEEGSAMIVQTSSGEYTIMSWLDRDVAENAARVSEGVVLAKYMPSEVLSLAVQTLYPKD"
        "GEIYLELSGYHREDLNFRKRRHHKRRARLRRASPSLQKSKELEEEKRLPLVPLDPEGAAA"
        "NLGSRRFRQRAHGAHHTHSRLLSRLSPSRYTSAGSSSVLQYSMSAECLVRAPDTFLFSFS"
        "INEASASAEPLNVLSGLPQRSSGLQAGSTEEEVARVWAKARIFKFPEMHCLFLHFNSLVF"
        "GTVTPILSLTITVLNLTHSSPLPRSIPTGIPLVKKTPSVLEKPKKVKAEPAKKKTNSQMP"
        "NYLEALNASSTEEVQALVSRLVLQPKPSRSVSPKKIDKAKAC"
    )

    # The UniProt sequence is 1-indexed. Python strings are 0-indexed.
    # TM3 ends at position 320. TM4 starts at position 432.
    # The linker is from position 321 to 431.
    # To slice this in Python, we use indices 320 to 431.
    linker_start_index = 320  # Corresponds to residue 321
    linker_end_index = 431   # Corresponds to the end of residue 431

    # Extract the TM3-TM4 linker sequence
    linker_sequence = full_sequence[linker_start_index:linker_end_index]

    # Count the number of Cysteine ('C') residues in the linker
    cysteine_count = linker_sequence.count('C')

    print("Human GABAA p1 Receptor (UniProt: P24046) Analysis")
    print("-" * 50)
    print("TM3-TM4 Linker Region: Residues 321-431")
    print("\nLinker Sequence:")
    # Print sequence wrapped for readability
    print(textwrap.fill(linker_sequence, 60))
    print("-" * 50)
    print(f"Number of Cysteine ('C') residues in the TM3-TM4 linker: {cysteine_count}")


if __name__ == '__main__':
    count_cysteines_in_linker()
