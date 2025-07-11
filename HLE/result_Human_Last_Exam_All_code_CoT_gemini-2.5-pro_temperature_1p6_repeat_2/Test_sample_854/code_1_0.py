def explain_lab_error():
    """
    Explains the critical error the laboratory made in trusting their Quality Control data for Batch 3 media.
    """
    print("The laboratory's mistake was a misinterpretation of their own Quality Control (QC) results. Here is the breakdown:")
    print("-" * 80)

    print("Step 1: Flawed Media Preparation in Batch 3")
    print("The antibiotic chloramphenicol was added to the PDA *before* it was autoclaved at 121 degrees for 25 minutes.")
    print("Chloramphenicol is heat-labile, meaning it is degraded and inactivated by the high heat of the autoclave process.")
    print("Result: Batch 3 was essentially PDA without any functional antibiotic.\n")

    print("Step 2: The Quality Control Test")
    print("The lab used Bacillus subtilis 6633, a bacterium, to test the quality of the agar.")
    print("The purpose of a QC test on media containing an antibiotic is to ensure the antibiotic *inhibits* bacterial growth.")
    print("Therefore, the correct 'expected result' should have been NO GROWTH of Bacillus subtilis on the plate.\n")

    print("Step 3: The Critical Error in Judgment")
    print("Because the chloramphenicol in Batch 3 was inactive, the Bacillus subtilis used for the QC test GREW on the plate.")
    print("The laboratory mistakenly interpreted this growth as an 'expected result', believing it meant the agar was fine.")
    print("This is the core mistake: They saw growth, which was evidence of a FAILED batch, but they interpreted it as a PASS.\n")

    print("Conclusion:")
    print("The laboratory believed the evidence because they fundamentally misunderstood what the evidence meant. The growth of bacteria on their QC plate for Batch 3 should have immediately signaled that the antibiotic was not working. Instead, they took it as a positive sign and proceeded to use the contaminated, non-selective media for their experiment.")

# Run the explanation
explain_lab_error()