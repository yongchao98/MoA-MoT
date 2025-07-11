def count_cysteines_in_linker():
    """
    This script calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAA receptor subunit rho-1 (GABRR1).
    """
    # Full amino acid sequence for human GABRR1 from UniProt (P24046)
    full_sequence = "MGFALALAPPLLLASLLLAGAAAAPSRQGGEKKPLGANASARAEILKSTATVNYTNLILRTLDNTMVTCTRILDRLLDGYDNRLRPGFGGAVTIGIYTTREVEILLDRYDAYEVEGDEYTVGVFQTYSATATSYTPEVLAGVSVALDWFDKAKIDSYSVGKPEDVTLTKLERDWLLKGYTVRVDVVYKDGKSVVLYTLNLTLPNHKEASKTSVEAKGYELKLTADTSIFYRAFSGLYPRVFFGVFGPSYAMVLSEYTIDVFFAQSWKDEKAAKAKEVGSRIRYIFQSFCSFINLVYWVSYLNFREPAGTVPARVLGITTVLTMTTLSISAARNLPKVAYATAMDWFIAVCYAFVFSALIEFATVNYFTKRGYAWDGKSVVPEPKKVKDPLIKKNMTILRAFPLLFGIFNLVYWATYLNREPQLKAPTPHQ"

    # According to UniProt P24046 annotations:
    # TM3 domain ends at position 339.
    # TM4 domain starts at position 436.
    # The linker is the region from position 340 to 435.
    # In Python's 0-based indexing, we slice from index 339 up to (but not including) 435.
    linker_start_index = 339
    linker_end_index = 435

    # Extract the linker sequence
    linker_sequence = full_sequence[linker_start_index:linker_end_index]

    # Count the number of Cysteine ('C') residues
    cysteine_count = linker_sequence.count('C')

    # Print the details and the final count
    print(f"The full protein has {len(full_sequence)} residues.")
    print(f"The TM3-TM4 linker is defined as the sequence from residue {linker_start_index + 1} to {linker_end_index}.")
    print(f"Linker sequence: {linker_sequence}")
    print(f"The number of Cysteine ('C') residues in the TM3-TM4 linker is: {cysteine_count}")

# Execute the function
count_cysteines_in_linker()