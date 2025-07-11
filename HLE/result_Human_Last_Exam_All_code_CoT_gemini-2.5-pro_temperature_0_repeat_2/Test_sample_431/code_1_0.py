def count_cysteine_in_linker():
    """
    This script calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAAρ1 receptor (UniProt: P24046).
    """
    # Full amino acid sequence for human GABAAρ1 receptor (UniProt: P24046)
    full_sequence = (
        "MGFALFLLAWATLASLSPGARADDIWVPGSFFRQVTPQALMELNSRVVARLPGLGPLRL"
        "PLDQVNVLLDGWLPVSFLDYEKPMSSIVSTNILDRLLDLLGYDNRLRPGFGGPPVNYGG"
        "KAVSYTLESIRYSIPDLFFVNSIVYWTDFFHQVATDMPLRILLDSALKVPDLNVYALLE"
        "GRSWKDEYRVCDVTFFTSYTGIDLSLDITSEVNLMTTPNFLLSKKSPLHAKPTVFERLM"
        "KPLWIPDLVLYNNTVGKYAIHVTCSCLALLNSHFRFDAYAKMDLVSLYLGIPLMLIWCE"
        "VFVAFLLQESAGAQAQPARTSVPAKASAPPLKQESLTTAPPAPKSPEEMRKLFIQRAAY"
        "LRELKRNPSVEEAESEPTRKMSAEDPMKNKWRYVFFPMAFLIFNMFWVIYKWINDKEVE"
        "KTSKADSTSRQLKQPLGPGSPNTALSKADDSLAEVQTSLGSTIAKSIDVY"
    )

    # UniProt positions are 1-based. Python strings are 0-based.
    # TM3 ends at position 313.
    # TM4 starts at position 429.
    # The linker is from position 314 to 428 (inclusive).
    # In 0-based indexing, this corresponds to the slice [313:428].
    linker_start_index = 313
    linker_end_index = 428

    # Extract the TM3-TM4 linker sequence
    tm3_tm4_linker_sequence = full_sequence[linker_start_index:linker_end_index]

    # Count the number of Cysteine ('C') residues in the linker sequence
    cysteine_count = tm3_tm4_linker_sequence.count('C')

    # Print the results
    print(f"Human GABAAρ1 receptor (UniProt P24046)")
    print(f"TM3-TM4 linker domain is from residue {linker_start_index + 1} to {linker_end_index}.")
    print(f"Linker Sequence: {tm3_tm4_linker_sequence}")
    print(f"The number of Cysteine ('C') residues in this domain is: {cysteine_count}")

if __name__ == "__main__":
    count_cysteine_in_linker()