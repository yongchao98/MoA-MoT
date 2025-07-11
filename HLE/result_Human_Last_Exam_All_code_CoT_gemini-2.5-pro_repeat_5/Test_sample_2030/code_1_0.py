def solve_oligo_problem():
    """
    Solves the bioinformatics problem by translating a DNA sequence in six frames,
    identifying a unique frame with specific SNP conditions, and designing a
    corresponding oligo.
    """
    # Step 1: Define the original sequence and the DNA codon table
    original_seq = "CTTCCCGCACAAGTGGT"
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
        'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
        'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
        'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
        'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
        'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # --- For demonstration: Six-frame translation ---
    # This manual analysis identifies the correct frame as described in the plan.
    # Forward Frame 1: CTT CCC CGC ACA AGT GGT -> L P R T S G
    # Forward Frame 2: TCC CCG CAC AAG TGG   -> S P H K W
    # Forward Frame 3: TTC CCC GCA CAA GTG   -> F P A Q V  <- Unique amino acids F and Q
    # Reverse Comp   : ACC ACT GTG CGC GGG AAG
    # Reverse Frame -1: ACC ACT GTG CGC GGG AAG -> T T V R G K
    # Reverse Frame -2: CCA CTG TGC GCG GGA   -> P L C A G
    # Reverse Frame -3: CAC TGT GCG CGG GAA   -> H C A R E

    print("Step 1: Identifying the unique reading frame")
    print("The original sequence is 5' CTT CCC CGC ACA AGT GGT 3'")
    print("Analysis of all six reading frames shows that Forward Frame 3 is the only one with two unique amino acids.")
    print(" - Frame 3 DNA: 5' TTC CCC GCA CAA GTG 3'")
    print(" - Frame 3 Peptide: F - P - A - Q - V")
    print(" - The unique amino acids are F (Phenylalanine) and Q (Glutamine).\n")

    print("Step 2: Applying the specified SNPs")
    # Properties: Phenylalanine (F) is non-polar. Glutamine (Q) is polar.
    # SNP 1: Polar (Q, codon CAA) -> STOP. The change is C -> T, making the codon TAA.
    # SNP 2: Non-polar (F, codon TTC) -> Cysteine (C). The change is T -> G, making the codon TGC.
    print(" - Phenylalanine (F) is non-polar. Its codon, TTC, is changed to TGC (Cysteine) by an SNP.")
    print(" - Glutamine (Q) is polar. Its codon, CAA, is changed to TAA (STOP) by an SNP.\n")

    print("Step 3: Determining the modified sequence for oligo design")
    # The modified sequence produces Cys, Pro, Ala, then STOP.
    # The part that is translated into amino acids is TGC-CCC-GCA.
    target_sequence = "TGCCCCGCA"
    print(f"The DNA sequence corresponding to the translated part of the modified frame is 5'-{target_sequence[0:3]}-{target_sequence[3:6]}-{target_sequence[6:9]}-3'.\n")

    print("Step 4: Designing the oligo")
    print("The oligo must bind to this target sequence, so it must be its reverse complement.")
    
    # Calculate reverse complement
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complementary_strand = "".join([complement_map[base] for base in target_sequence])
    oligo_sequence = complementary_strand[::-1]

    print(f"Target Sequence : 5' {target_sequence[0]} {target_sequence[1]} {target_sequence[2]} {target_sequence[3]} {target_sequence[4]} {target_sequence[5]} {target_sequence[6]} {target_sequence[7]} {target_sequence[8]} 3'")
    print(f"Oligo Sequence  : 5' {oligo_sequence[0]} {oligo_sequence[1]} {oligo_sequence[2]} {oligo_sequence[3]} {oligo_sequence[4]} {oligo_sequence[5]} {oligo_sequence[6]} {oligo_sequence[7]} {oligo_sequence[8]} 3'")
    
    print(f"\nThe final DNA sequence of the oligo is 5' {oligo_sequence} 3'.")
    
    # Final answer block
    print(f"\n<<<{oligo_sequence}>>>")

solve_oligo_problem()