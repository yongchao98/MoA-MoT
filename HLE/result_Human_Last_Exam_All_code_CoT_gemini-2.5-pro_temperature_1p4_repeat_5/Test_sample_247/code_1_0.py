import textwrap

def explain_tcr_sequencing_plan():
    """
    This function explains the necessary steps to modify a 3' single-cell
    RNA-seq workflow to successfully capture and sequence TCR CDR3 regions.
    """
    
    # Header
    print("="*60)
    print("Plan for Single-Cell TCR CDR3 Sequencing".center(60))
    print("="*60)

    # The Problem
    problem_text = """The existing method uses a poly(dT) primer to capture the 3' end of mRNA. The target CDR3 region is at the 5' end of the TCR transcript, which is too far to be captured and sequenced reliably with this 3'-biased approach."""
    print("\n[Problem Analysis]")
    print(textwrap.fill(problem_text, width=60))

    # The Solution
    solution_text = """The solution is to introduce a targeted amplification step after cDNA synthesis. This involves designing specific PCR primers that flank the CDR3 region to selectively enrich for the TCR sequences of interest, creating a dedicated TCR library for sequencing."""
    print("\n[Proposed Solution]")
    print(textwrap.fill(solution_text, width=60))

    # Detailed Steps
    print("\n[Actionable Steps]")
    
    # Step 1: cDNA Synthesis
    print("\n1. First, ensure cDNA containing the CDR3 region is synthesized.")
    step1_details = """While the original poly(dT) capture is inefficient, a robust protocol would involve priming the reverse transcription from the TCR constant (C) region. This can be done by adding TCR-specific primers to the reverse transcription reaction."""
    print(textwrap.fill(step1_details, width=58, initial_indent="  ", subsequent_indent="  "))
    
    # Step 2: Targeted PCR Amplification
    print("\n2. Second, perform targeted PCR to amplify the V(D)J region.")
    step2_details = """This is the most critical step. From the total cDNA pool, you must specifically amplify the TCR transcripts. Design PCR primers that 'flank' the CDR3 region:"""
    print(textwrap.fill(step2_details, width=58, initial_indent="  ", subsequent_indent="  "))
    
    # Primer design details
    print(textwrap.fill("- Forward Primers: A cocktail of primers targeting the conserved portions of the different V-gene families (upstream of CDR3).", width=56, initial_indent="    ", subsequent_indent="    "))
    print(textwrap.fill("- Reverse Primer(s): One or more primers targeting the TCR constant (C) region (downstream of CDR3).", width=56, initial_indent="    ", subsequent_indent="    "))
    
    # Explanation of PCR outcome
    pcr_outcome_text = """This PCR strategy specifically enriches for the desired fragments, making them the vast majority of the final sequencing library. A nested PCR approach (using two rounds of PCR with inner and outer primers) is often used to increase specificity."""
    print(textwrap.fill(pcr_outcome_text, width=58, initial_indent="  ", subsequent_indent="  "))
    
    # Conclusion
    print("\n[Final Outcome]")
    final_text = """This modified protocol will produce a library highly concentrated with TCR CDR3 sequences, which can then be successfully sequenced using the existing 75x225 Illumina modality."""
    print(textwrap.fill(final_text, width=60))
    print("\nThis strategy directly corresponds to designing primers that bind to regions flanking the CDR3 for subsequent PCR amplification.")


# Execute the function to print the explanation.
explain_tcr_sequencing_plan()