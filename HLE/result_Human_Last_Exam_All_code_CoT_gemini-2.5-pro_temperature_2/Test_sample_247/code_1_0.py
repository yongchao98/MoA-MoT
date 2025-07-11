import textwrap

def simulate_tcr_sequencing():
    """
    Simulates TCR sequencing to demonstrate why 3' capture fails for CDR3
    and how targeted PCR amplification (Option D) solves the problem.
    """
    # 1. Define the structure of a model TCR-beta transcript (lengths are approximate)
    # The CDR3 region is part of the V-D-J junction.
    regions = {
        "5UTR": "-" * 100,
        "Leader": "L" * 50,
        "V_Region": "V" * 300,
        "CDR3_Region": "CDR3" * 10,  # "CDR3CDR3..." total 40 chars
        "J_Region": "J" * 50,
        "C_Region": "C" * 400,
        "3UTR": "=" * 200,
        "PolyA_Tail": "A" * 30
    }
    tcr_mrna = "".join(regions.values())
    mrna_len = len(tcr_mrna) - len(regions["PolyA_Tail"])
    
    # Define bead oligo and read lengths from the problem
    bead_oligo = "BARCODE_UMI_" * 5 # Represents the 50bp oligo with barcode and UMI
    read1_len = 75
    read2_len = 225

    print("--- TCR mRNA & System Parameters ---")
    print(f"Model TCR mRNA length (excluding PolyA tail): {mrna_len} bases")
    print(f"Sequencing: Read 1 = {read1_len}bp, Read 2 = {read2_len}bp\n")
    print("Full mRNA structure:")
    print(textwrap.fill(f"5' -> {tcr_mrna} -> 3'", width=80), "\n")
    
    # 2. Simulate the PhD student's original method (3' Gene Expression)
    print("--- Scenario 1: Original 3' Capture Method ---")
    print("Process: Reverse transcription is primed from the PolyA tail.")
    
    # In 3' sequencing, Read 2 sequences the transcript from the 3' end.
    read2_content_original = tcr_mrna[-(read2_len + len(regions["PolyA_Tail"])):-len(regions["PolyA_Tail"])]
    
    print(f"Read 2 (225bp) sequences the following part of the transcript:")
    print(textwrap.fill(f"...{read2_content_original}", width=80))
    
    # Check if CDR3 is in the sequenced part
    if regions["CDR3_Region"] in read2_content_original:
        print("\nResult: SUCCESS - CDR3 sequence was captured. (This should not happen)")
    else:
        print("\nResult: FAILURE - The 225bp read from the 3' end is too short to reach the CDR3 region.")
        print("The read only covers the 3'UTR and a portion of the Constant (C) region.\n")

    # 3. Simulate the proposed solution from Option D (Targeted PCR)
    print("--- Scenario 2: Proposed Solution (Option D) ---")
    print("Process: 1. Create full-length cDNA library (same as before).")
    print("         2. Use PCR to specifically amplify the V(D)J region.")
    
    # The full cDNA is synthesized from the polyA tail, attached to the bead oligo
    # For simplicity, we model the PCR product directly.
    # Forward primer binds to bead oligo sequence.
    # Reverse primer binds within the Constant (C) region.
    # This creates a smaller amplicon containing the 5' end of the transcript.
    
    pcr_product_end = (len(regions["5UTR"]) + len(regions["Leader"]) + 
                       len(regions["V_Region"]) + len(regions["CDR3_Region"]) + 
                       len(regions["J_Region"]) + len(regions["C_Region"]))
    
    pcr_amplicon = bead_oligo + tcr_mrna[:pcr_product_end]

    print("\nA targeted PCR amplicon is created. It looks like this:")
    print(textwrap.fill(f"5' -> {pcr_amplicon[:60]}...{pcr_amplicon[-60:]} -> 3'", width=80))
    print(f"Total Amplicon Length: {len(pcr_amplicon)} bp\n")

    print("Now, we sequence this new, shorter PCR product:")
    # Read 1 still sequences the barcode/UMI from the bead oligo end.
    read1_content_new = pcr_amplicon[:read1_len]
    print("Read 1 (75bp) captures the barcode and UMI from the bead oligo part.")
    
    # Read 2 sequences from the other end of the amplicon, which is the V(D)J region.
    # Note: In reality, sequencing adapters are ligated, and Read 2 starts from a sequencing primer
    # binding site on the adapter. It then reads into the insert (V region first).
    
    insert_to_sequence_for_read2 = tcr_mrna[:pcr_product_end]
    read2_content_new = insert_to_sequence_for_read2[:read2_len]
    
    print(f"\nRead 2 (225bp) sequences the other end of the insert:")
    print(textwrap.fill(f"{read2_content_new}...", width=80))
    
    # Check if CDR3 is in the new Read 2
    if regions["CDR3_Region"] in read2_content_new:
        print(f"\nResult: SUCCESS! The CDR3 region ('{regions['CDR3_Region']}') is present in the sequence from Read 2.")
    else:
        print("\nResult: FAILURE - CDR3 sequence was still not captured.")

if __name__ == "__main__":
    simulate_tcr_sequencing()