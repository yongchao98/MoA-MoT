def solve_cellular_automata_puzzle():
    """
    This function holds the derived solution for the cellular automata puzzle.
    The logic for deriving the groups is explained in the text above.
    """
    R1_Labels = "GCKE"
    R2_Labels = "FHMA"
    R3_Labels = "DNLP"
    R4_Labels = "IBOJ"

    # Presenting the answer in the required format {R1_Labels, R2_Labels, R3_Labels, R4_Labels}
    final_answer = f"{{{R1_Labels},{R2_Labels},{R3_Labels},{R4_Labels}}}"
    print(final_answer)

solve_cellular_automata_puzzle()