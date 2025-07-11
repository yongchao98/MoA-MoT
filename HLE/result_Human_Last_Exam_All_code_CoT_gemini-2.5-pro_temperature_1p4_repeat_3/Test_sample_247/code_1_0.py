def design_tcr_enrichment_primers():
    """
    This function demonstrates the primer design strategy for enriching TCR
    sequences from a single-cell cDNA library.
    """

    # The student's BD Rhapsody beads have a universal PCR handle. The sequence
    # for the Rhapsody system's library amplification is publicly known.
    # This primer will bind to the universal handle on the bead-oligo side of the cDNA.
    pcr_handle_sequence = "AAGCAGTGGTATCAACGCAGAGT"
    reverse_primer = pcr_handle_sequence

    # For the forward primer, we need a sequence from the constant region of the
    # T-cell receptor transcript. This region is 3' of the V(D)J (CDR3) region,
    # so a primer here will ensure the entire CDR3 is copied during PCR.
    # We will use a conserved part of the human TCR Beta Constant Region 1 (TRBC1).
    tcr_beta_constant_region_sequence = "GACCTCAGATCCAGTTTAAGGATTCCTTCACTGATATCCTGCAGATCACAGAGTCCTCTGTCTACCAGTTTGATTCT"
    
    # We can select a primer from the beginning of this sequence.
    # A typical PCR primer is ~20-25 nucleotides long.
    forward_primer = tcr_beta_constant_region_sequence[0:25]

    print("To solve the issue, targeted PCR enrichment should be performed.")
    print("This requires two specific primers for a PCR reaction following reverse transcription:\n")
    
    print("1. Forward Primer (TCR Constant Region specific):")
    print("   - Target: TCR Beta Constant Region (TRBC)")
    print(f"   - Sequence (5' to 3'): {forward_primer}")
    print("\n")
    
    print("2. Reverse Primer (Bead Oligo specific):")
    print("   - Target: Universal PCR handle on the bead oligo")
    print(f"   - Sequence (5' to 3'): {reverse_primer}")
    print("\n")

    print("--- PCR Strategy Explanation ---")
    print("The cDNA from the initial reaction serves as a template.")
    print("The Forward Primer binds to the TCR constant region on the cDNA.")
    print("The Reverse Primer binds to the universal PCR handle at the other end of the cDNA.")
    print("PCR amplification will selectively create many copies of the full TCR sequence,")
    print("including the V(D)J and CDR3 regions, making them easy to find during sequencing.")
    
# Execute the function to display the primer design strategy.
design_tcr_enrichment_primers()
