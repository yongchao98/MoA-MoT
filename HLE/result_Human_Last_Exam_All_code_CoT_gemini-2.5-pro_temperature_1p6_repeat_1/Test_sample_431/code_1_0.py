import sys

def count_cysteine_in_linker():
    """
    This function calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAA receptor subunit rho-1 (GABRR1).

    The sequence and domain information are based on the UniProt entry P24046.
    """
    # Canonical sequence for human GABRR1 from UniProt (P24046)
    protein_sequence = (
        "MGFHAGRRPSLPSLPLLLLLLLLLRAGPAQGSEVYVGSYQVDLNDRGYDNRLRPGFGGAV"
        "TIEVTDIFAYTNIFLFQQWRDNRLKYNTNVLVQTMYNNIMSKSTVPLNMLDIPMDAMFYT"
        "PKAYAPCASFDIRYGYCTDFEFFKGSAAPRFTFGIVNSVLLSLSFEVSEVDSGSRLNQYE"
        "DLLGQTVRIGNYIVWTGFAPYLVFEMSVGKSPIEGTTMDYFRYNAGYKPCLGVCYALAFW"
        "YAFVANYLPLDMSYEVSKKTVEFTRSGSYAYITEDMFFVTCFFLRLALDWPYAYAIKCWI"
        "PDLFFANSKSADPHYTPSTMSKAPQSPSAKKTEEPPESAPASKKTEEPPQSKPTAPPEPA"
        "KKTVPAPAVSKEEAPSPPKKAEEAEPAKTEAPPAKKAEEASPAKKAEDAVPTKKAEAAAQ"
        "KTAAAGPLHQPNGPAGKIPQNPDQKRLEGSRIGFWLPLVLYFLAANMMVAEVAYCLHDHG"
        "ALASPPSQLKAPKSLLPGRSVKIPIPYALCVFNVLWIIWVSFVSRFPIPF"
    )

    # According to UniProt annotation for P24046:
    # Transmembrane domain 3 (TM3) ends at position 335.
    # Transmembrane domain 4 (TM4) starts at position 436.
    tm3_end_pos = 335
    tm4_start_pos = 436

    # The linker is the region between TM3 and TM4.
    # We use 1-based indexing for residue positions.
    linker_start_pos = tm3_end_pos + 1
    linker_end_pos = tm4_start_pos - 1

    # In Python, string indexing is 0-based, and slicing is exclusive of the end index.
    # So, to get residues from linker_start_pos to linker_end_pos (inclusive),
    # we slice from (linker_start_pos - 1) to linker_end_pos.
    linker_sequence = protein_sequence[linker_start_pos - 1 : linker_end_pos]

    # Count the number of Cysteine ('C') residues in the linker sequence
    cysteine_count = linker_sequence.count('C')

    # Output the results
    print(f"Protein: Human GABAA rho-1 receptor (UniProt: P24046)")
    print(f"TM3 domain ends at position: {tm3_end_pos}")
    print(f"TM4 domain starts at position: {tm4_start_pos}")
    print(f"The TM3-TM4 linker is defined as the region from residue {linker_start_pos} to {linker_end_pos}.")
    print(f"Linker sequence: {linker_sequence}")
    print(f"\nFinal Equation:")
    print(f"Number of Cysteine residues in TM3-TM4 linker (residues {linker_start_pos}-{linker_end_pos}) = {cysteine_count}")

# Execute the function
count_cysteine_in_linker()

# We use sys.stdout.flush() to ensure the output is printed before the marker.
sys.stdout.flush()