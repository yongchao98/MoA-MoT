import textwrap

def demonstrate_tcr_amplification():
    """
    This script demonstrates the solution (Option D) for capturing the TCR CDR3 region.
    The strategy is to use targeted PCR after initial cDNA synthesis. A forward primer
    in the TCR constant region and a reverse primer on the bead adapter sequence are used
    to amplify a specific fragment containing the V(D)J/CDR3 region.
    """
    print("--- Step 1: Define the components of the TCR mRNA and the bead adapter ---")

    # Simplified representation of the universal adapter on the sequencing bead.
    # In a real experiment, this region would also contain the cell label and UMI.
    ADAPTER = "GATTACAGATTACAGATTACAGATTACA"

    # The components of a sample TCR Beta chain mRNA, from 5' to 3' end.
    TCR_MRNA_PARTS = {
        "5_UTR_and_Leader": "AGCTGTGGAGTCC...",
        "V_GENE": "...ATGTCTGAGAGACCCAGAG...",
        "CDR3": "TGTGCCAGCAGCTTAG...", # The critical region of interest
        "J_GENE": "...ACAGGGGCAGGGACACAG...",
        "C_GENE": "GAGCAGTTCTTCGGGCCAGGGACACGGCTCACCTGT...", # Constant region
        "3_UTR_and_PolyA": "...GGGTCCTGAGTGACACAAAAAAAAAAAA"
    }
    
    # Print the components for clarity
    print(f"Bead Adapter Sequence: {ADAPTER}")
    print("TCR mRNA Components:")
    for name, seq in TCR_MRNA_PARTS.items():
        print(f"  - {name}: {seq}")
    print("-" * 50)


    print("\n--- Step 2: Design primers for targeted PCR amplification (The Solution) ---")
    
    # The forward primer binds within the constant region (C_GENE).
    # This primes amplification towards the 5' end of the transcript, through the CDR3.
    FWD_PRIMER = TCR_MRNA_PARTS["C_GENE"][0:20]

    # The reverse primer is complementary to the adapter sequence from the bead.
    # In practice, a primer with the same sequence as the adapter is used to bind
    # the reverse-complement strand of the cDNA.
    REV_PRIMER = ADAPTER

    print(f"Forward Primer (from Constant Region): {FWD_PRIMER}")
    print(f"Reverse Primer (from Bead Adapter): {REV_PRIMER}")
    print("-" * 50)


    print("\n--- Step 3: Construct the final PCR product (Amplicon) ---")
    print("This amplicon is the result of the targeted PCR. It is a short DNA fragment")
    print("flanked by the primer sequences, containing the complete V(D)J region.")

    # The final amplicon is constructed from the relevant parts of the mRNA,
    # flanked by the primer sequences.
    amplicon_components = [
        FWD_PRIMER,  # Start of the amplicon
        "..." + TCR_MRNA_PARTS["J_GENE"],
        TCR_MRNA_PARTS["CDR3"],
        TCR_MRNA_PARTS["V_GENE"],
        TCR_MRNA_PARTS["5_UTR_and_Leader"] + "...",
        REV_PRIMER # End of the amplicon
    ]
    final_amplicon = "".join(amplicon_components)

    print("\nFinal Amplicon Structure:")
    print(f"  [Forward Primer]...[J]...[CDR3]...[V]...[5' UTR]...[Reverse Primer]")
    
    print("\nFinal Amplicon Sequence (wrapped for readability):")
    print("\n".join(textwrap.wrap(final_amplicon, width=80)))
    print("-" * 50)


    print("\n--- Conclusion ---")
    print("The CDR3 region is now successfully captured within a relatively short amplicon.")
    print("With a 75x225 sequencing run:")
    print(" - Read 1 (75bp) sequences the 'Reverse Primer' end, capturing the Adapter, Cell Label, and UMI.")
    print(" - Read 2 (225bp) sequences from the 'Forward Primer' end, easily covering the C, J, CDR3, and V regions.")
    print("This solves the student's problem without changing the sequencer or beads.")

if __name__ == '__main__':
    demonstrate_tcr_amplification()