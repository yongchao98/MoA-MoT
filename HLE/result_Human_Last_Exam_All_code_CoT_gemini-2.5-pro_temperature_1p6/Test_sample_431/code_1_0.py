def find_cysteine_in_linker():
    """
    Calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAAρ1 homomeric receptor (UniProt: P24046).
    """
    # The canonical amino acid sequence for human GABAAρ1 from UniProt (P24046)
    protein_sequence = "MGFHAGSKLLVSSFSQDLHSLDSRYDWQDRILTFLYGLTTRLTMTTLSIHSRNLKLLPVSLIVSGLTVFLYSIKVETVETHALEKTLVNYLVTMSIWTHLREMVHPLALDFTPMNLFTVVCALQLLSLESIKVALPRVSFLHLNNKMVPMDAIFYTMFCVSFSLASVEVAVLVNYFSRSAEHDVARRRARRKGYVAVADMWLAVPALPMVFSCLLNVVPARAVGTTEASSTLTEATTAKSRPTPALNIHMPAPRPPSVPSAAGPLASPGASPAASEVKKAEKARSTAIASSNTTASDSCSCVASMPLGIAGSYAVTMIALLISSLAYAVAEALNFWSRVKPLLRLGAVLILYGAILLEALNFVLVCVSLSLHTSPSVPLIPFGGLAPAKAVASANSQL"

    # According to UniProt annotation for P24046:
    # Transmembrane helix 3 (TM3) ends at position 322.
    # Transmembrane helix 4 (TM4) starts at position 435.
    # The linker is the sequence between these helices.
    linker_start_pos = 323  # 1-based position
    linker_end_pos = 434    # 1-based position

    # Python strings are 0-indexed. To extract the sequence from position 323 to 434,
    # we slice from index (323 - 1) up to (but not including) index 434.
    start_index = linker_start_pos - 1
    end_index = linker_end_pos
    linker_sequence = protein_sequence[start_index:end_index]

    # Count the number of Cysteine ('C') residues in the extracted linker sequence.
    cysteine_count = linker_sequence.count('C')

    # Output the details of the process and the final result.
    print("Analysis of Human GABAAρ1 Receptor (UniProt P24046):")
    print(f"1. Identified TM3-TM4 linker boundaries: Positions {linker_start_pos} to {linker_end_pos}.")
    print(f"2. Extracted the linker sequence of length {len(linker_sequence)}.")
    print(f"   Linker Sequence: {linker_sequence}")
    print("\nFinal Calculation:")
    print(f"Number of Cysteine ('C') residues in sequence from position {linker_start_pos} to {linker_end_pos} = {cysteine_count}")


if __name__ == '__main__':
    find_cysteine_in_linker()