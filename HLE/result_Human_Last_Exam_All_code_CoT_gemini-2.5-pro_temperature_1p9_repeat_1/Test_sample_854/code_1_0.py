import textwrap

def explain_lab_error():
    """
    This function prints a step-by-step explanation for the mistake made by the laboratory.
    """
    
    # Define the explanation steps
    title = "Analysis of the Laboratory's Mistake"
    line = "=" * len(title)

    step1_title = "1. The Critical Difference in Preparation:"
    step1_text = """The preparation of Batch 3 differed from Batches 1 and 2 in one critical step. Person B added the antibiotic chloramphenicol to the media *before* the autoclaving process. The procedure states: 'the PDA batch was autoclaved at 121 degrees under high pressure for 25 minutes' after the chloramphenicol was added."""
    
    step2_title = "2. The Chemical Property of the Antibiotic:"
    step2_text = """Chloramphenicol is a heat-labile antibiotic. This means it is sensitive to heat and breaks down or becomes inactive when subjected to high temperatures like those used in an autoclave (121 degrees Celsius). The standard and correct procedure is to autoclave the media first, let it cool to a safe temperature (around 50°C), and then add any heat-sensitive supplements like chloramphenicol aseptically."""

    step3_title = "3. The Consequence of the Error:"
    step3_text = """By autoclaving the media after adding the antibiotic, the chloramphenicol in Batch 3 was destroyed. This left the agar with no effective agent to prevent bacterial growth. Batches 1 and 2, which were presumably prepared correctly with the antibiotic added after autoclaving, retained their antibacterial properties."""

    step4_title = "4. Why the QC Check Failed and Evidence Was Misleading:"
    step4_text = """The lab believed Batch 3 was safe because it was prepared in a biosafety cabinet, which seemed more rigorous, and it passed a Quality Control (QC) check. The QC check was flawed; it failed to detect the absence of antibacterial activity. This could be because the QC organism (*Bacillus subtilis*) was in its spore form and didn't germinate during the test, or the lab simply misinterpreted the results."""
    
    conclusion_title = "5. Final Conclusion:"
    conclusion_text = """The laboratory's mistake was believing their evidence because they trusted a seemingly superior preparation process and a flawed QC check. They overlooked a fundamental principle of microbiology: the heat-labile nature of chloramphenicol. When all batches were exposed to common airborne bacteria (like the spore-forming rods found), only Batch 3 lacked the protection to inhibit bacterial growth, leading to the observed contamination."""

    # Print the formatted explanation
    print(f"{title}\n{line}\n")
    
    print(step1_title)
    print(textwrap.fill(step1_text, width=80))
    print("\n" + "-"*80 + "\n")

    print(step2_title)
    print(textwrap.fill(step2_text, width=80))
    print("\n" + "-"*80 + "\n")

    print(step3_title)
    print(textwrap.fill(step3_text, width=80))
    print("\n" + "-"*80 + "\n")

    print(step4_title)
    print(textwrap.fill(step4_text, width=80))
    print("\n" + "-"*80 + "\n")
    
    print(conclusion_title)
    print(textwrap.fill(conclusion_text, width=80))

# Execute the function to print the analysis
explain_lab_error()
