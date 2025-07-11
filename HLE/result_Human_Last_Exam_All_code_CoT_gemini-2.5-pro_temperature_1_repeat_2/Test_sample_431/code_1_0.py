def count_cysteine_in_linker():
    """
    This function calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAA rho-1 receptor (UniProt: P24046).
    """
    # Full protein sequence for human GABAA rho-1 receptor (GABRR1_HUMAN, P24046)
    full_sequence = (
        "MGFALWAVWALWVLGSRVSEVQDQLPDSSTINYTMNDLFLDWNPADATSYTPNLSLAVPLD"
        "GSVTVEMTHMTLVKSTYIPDTYFLNDKKSFVHGVTVKNRMIRLHPDGTVLYGLRITTTAAC"
        "MMDLRRYPLDEQNCTLEIESYGYTTDDIEFYWRGGDKAVTGVERIELPQFSIVEHRLVSRR"
        "FGYATLNLSFVPFSDSYVNNDMVSFCLFDPMYAFAMAFNVLVYEWLIGRPLTVPATVGVTT"
        "KLSLPIPSDLGINAILSYIASCTTIMGVWINSGHIPLAKSPIKRALPPNGKSPASIKTTET"
        "IKAYAKAEVKKKSAAKKIVKDAHIKKLEDRASKEMARELIKEVMTRFASLFNHEEGETTGE"
        "GALVPGYFLIERKSLKIKTIGYFVIFFSLSFYLATTLYYPLKIPSTVIVSWVSFWINMDAA"
        "PARVALGITTVLTMTTISTHIALQRLPFHLFHFKRNIGYFILQTYLPCIMTVILSQVSFWL"
        "NRESVPARTVFGVTTVLTMTTLSISARNSLPKVAYATAMDWFIAVCYAFVFSALIEFATVN"
        "YFTKRAE"
    )

    # According to UniProt P24046, the TM3-TM4 linker (a large cytoplasmic domain)
    # is located from amino acid position 322 to 433.
    # Note: Protein sequence positions are 1-based, while Python string indices are 0-based.
    # To get the sequence from 322 to 433, we slice from index 321 to 433.
    linker_start_index = 322 - 1
    linker_end_index = 433

    # Extract the linker sequence
    tm3_tm4_linker_sequence = full_sequence[linker_start_index:linker_end_index]

    # Count the number of Cysteine ('C') residues in the linker sequence
    cysteine_count = tm3_tm4_linker_sequence.count('C')

    # Print the result
    print(f"The sequence of the TM3-TM4 linker is: {tm3_tm4_linker_sequence}")
    print(f"The number of Cysteine residues found in the TM3-TM4 linker is: {cysteine_count}")

# Execute the function
count_cysteine_in_linker()