def count_cysteine_in_linker():
    """
    This script finds the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAArho1 receptor (NCBI RefSeq NP_002035.1).
    """

    # 1. Full amino acid sequence for human GABAA rho-1 receptor from NCBI NP_002035.1
    # Note: Sequence is joined from its FASTA representation.
    full_sequence = (
        "MSSSPPQPSQAKETNLDSNRTVLYSLSIGVSSVSETFHDAKVTPKGTPVPAAKESLDLLGYTLETIYQND"
        "TGVEYAVDVEFWRSVDYALPLEFGGYMKFIQRAWYDERYDVFFRQQWSDERLKFKGPVQRVPAGPNMLLM"
        "TSISIKASIDMVSEVNVMDYTLTMYFQQAWRDRRLAAYAVGMGILLALAYSFSALIEAATVNYFTKRGW"
        "AWDGKSVVNDAKGRMLRLIWQPDLFFVNFSSIPVLALFSVNLVYWATYLNFLRKTVPDAYILCTLWILDY"
        "QGSSALPSSLAALVALSYLTMWFIAWFLSRPEDVASTIIKSKGAAPRSLSAMVTMTLLSISSARNSLPKV"
        "TYATVMDWFIAVCYAFVFSALIEFATVNYFTKRGYAWDGNKKIGFPINFIAFVDTVEIAPSMGIPDGEAT"
        "AELPDTPCDARASGRPLLPQSSPGPLQDLQDPSLCASLLGEDSPLPHLPAFSG"
    )

    # 2. Define the linker boundaries based on NCBI annotation for NP_002035.1
    # TM3: 297-317
    # TM4: 429-449
    # Linker is from residue 318 to 428.
    linker_start_residue = 318
    linker_end_residue = 428

    # 3. Extract the linker sequence. Python uses 0-based indexing.
    # To get residue N, we use index N-1.
    # To include the end residue, the slice goes up to end_residue's index + 1.
    start_index = linker_start_residue - 1
    end_index = linker_end_residue

    linker_sequence = full_sequence[start_index:end_index]

    # 4. Find the positions of all Cysteine ('C') residues in the linker.
    cysteine_positions = []
    # Enumerate through the linker sequence to find Cysteines
    for i, amino_acid in enumerate(linker_sequence):
        if amino_acid == 'C':
            # The position in the full sequence is the starting index of the linker
            # plus the index within the linker, plus 1 (for 1-based residue numbering).
            position_in_full_sequence = start_index + i + 1
            cysteine_positions.append(position_in_full_sequence)
    
    cysteine_count = len(cysteine_positions)

    # 5. Print the results clearly.
    print(f"Analysis of Human GABAAρ1 Receptor (NP_002035.1)")
    print(f"---------------------------------------------------")
    print(f"The TM3-TM4 linker domain spans residues: {linker_start_residue} to {linker_end_residue}")
    print(f"Linker Sequence: {linker_sequence}")
    
    if cysteine_count > 0:
        # Build the "equation" part of the output string
        equation_str = " + ".join([f"C at pos {pos}" for pos in cysteine_positions])
        print(f"\nFound Cysteine residues at the following positions: {cysteine_positions}")
        print(f"Calculation: {equation_str}")
        print(f"\nFinal count of Cysteine residues in the linker: {cysteine_count}")
    else:
        print("\nNo Cysteine residues were found in the TM3-TM4 linker domain.")


if __name__ == "__main__":
    count_cysteine_in_linker()
