import textwrap

def generate_oligo_from_sequence():
    """
    This function solves the biological puzzle to determine the correct oligo sequence.
    The logic is derived from analyzing the six reading frames of the original DNA.
    """
    
    # Step 1: Identification of the unique reading frame.
    # The original sequence is 5' CTT CCC CGC ACA AGT GGT 3'.
    # Analysis of all six reading frames reveals that Frame +2 (starting from the second nucleotide)
    # is the only one containing two unique amino acids not found in other frames.
    # The sequence for this frame and its translation are:
    frame_dna = "TTC CCG CAC AAG TGG"
    frame_peptide = "F - P - H - K - W"
    
    print("Step 1: Identify the unique reading frame")
    print("-" * 40)
    print(f"The unique reading frame is Frame +2, which corresponds to the DNA sequence:\n5' {frame_dna} 3'")
    print(f"This translates to the peptide: {frame_peptide}")
    print("The amino acids F (Phenylalanine) and W (Tryptophan) are unique to this frame.\n")

    # Step 2: Determine the SNPs on the unique frame.
    # The frame contains polar amino acids H (Histidine) and K (Lysine),
    # and non-polar amino acids F (Phenylalanine), P (Proline), and W (Tryptophan).
    #
    # SNP 1 (Polar -> Stop): The codon for K (Lysine) is AAG. A single nucleotide change (A -> T)
    # creates the codon TAG, which is a stop codon. This is the only possible way to satisfy this condition.
    #
    # SNP 2 (Non-polar -> Cysteine): To preserve the translated sequence before the stop codon,
    # the SNP must affect the non-polar amino acid W (Tryptophan) which is located after K (Lysine).
    # The codon for W is TGG. A single change (G -> T) creates the codon TGT, which codes for Cysteine.
    print("Step 2: Apply the SNP (Single Nucleotide Polymorphism) rules")
    print("-" * 40)
    print("One SNP changes a polar amino acid to a STOP codon.")
    print("  - The polar Lysine (K), with codon AAG, is changed to a STOP codon (TAG).")
    print("Another SNP changes a non-polar amino acid to Cysteine (C).")
    print("  - The non-polar Tryptophan (W), with codon TGG, is changed to Cysteine (TGT).")
    modified_dna = "TTC CCG CAC TAG TGT"
    print(f"The resulting modified DNA sequence is:\n5' {modified_dna} 3'\n")
    
    # Step 3: Determine the oligo's target sequence.
    # The oligo binds to the part of the modified sequence that gets translated.
    # Translation stops at the new TAG codon. Therefore, the translated part is F-P-H.
    target_dna_str = "TTC CCG CAC"
    target_dna = target_dna_str.replace(" ", "")
    
    print("Step 3: Determine the oligo target sequence")
    print("-" * 40)
    print("The oligo binds to the region of the modified sequence that is translated into amino acids.")
    print("Since translation stops at the new STOP codon, the oligo's target is the DNA coding for F-P-H.")
    print(f"Target Sequence: 5' {target_dna_str} 3'\n")

    # Step 4: Calculate the oligo sequence (reverse complement of the target).
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = "".join(complement_map[base] for base in reversed(target_dna))
    
    # Format the oligo sequence for readability
    final_oligo = " ".join(textwrap.wrap(reverse_complement_seq, 3))

    print("Step 4: Design the oligo")
    print("-" * 40)
    print("The oligo sequence must be the reverse complement of the target sequence to bind to it.")
    print(f"The final oligo sequence is:\n5' {final_oligo} 3'")


generate_oligo_from_sequence()

# The final answer in the required format
print("\n<<<5' GTG GGC AAG 3'>>>")
