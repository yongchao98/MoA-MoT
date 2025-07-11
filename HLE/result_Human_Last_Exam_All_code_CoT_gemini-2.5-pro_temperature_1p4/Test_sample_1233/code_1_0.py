import textwrap

def solve_microbiology_case():
    """
    Analyzes a clinical microbiology scenario to determine the best corrective action.
    """
    # Define the core facts from the case study
    facts = {
        "Symptom": "Bloody diarrhea",
        "Suspected Pathogen": "Campylobacter jejuni (based on selective media, 42C, and microaerophilic conditions)",
        "Pathogen Characteristics": "Can be slow-growing, especially if initial sample viability is low.",
        "Lab 1 Action": "Incubated plates for 'two days' (48 hours).",
        "Lab 1 Result": "Bacillus species (a likely contaminant or overgrowth), a presumptive negative for Campylobacter."
    }

    # The central problem to solve
    problem = "The target organism, Campylobacter, was likely present but not detected. This can happen if the organism is slow to grow or present in low numbers, and is outcompeted by faster-growing contaminants."

    # The options provided
    choices = {
        'A': "Obtain a fresh sample",
        'B': "Decrease sample processing time",
        'C': "Used Fresh Campylobacter plates",
        'D': "Incubated the sample for longer",
        'E': "Increase the number of plates"
    }

    # Print the analysis step-by-step
    print("Analyzing the Laboratory Scenario:")
    print("-" * 35)
    print(f"The patient's symptom ({facts['Symptom']}) and the lab's specific incubation protocol point towards a suspected infection with {facts['Suspected Pathogen']}.")
    print(f"The first lab's result ({facts['Lab 1 Result']}) indicates they failed to isolate the true pathogen.")
    print(f"\nProblem Statement: {textwrap.fill(problem, 70)}")

    print("\nEvaluating the Potential Solutions:")
    print("-" * 35)

    # A
    print(f"A. {choices['A']}:")
    print("   This would likely lead to a correct diagnosis, but it bypasses the problem with the original sample rather than solving it.\n")

    # B
    print(f"B. {choices['B']}:")
    print("   This is an important preventative measure for future samples, as delays can decrease Campylobacter viability. However, it doesn't fix the plates that have already been inoculated.\n")

    # C
    print(f"C. {choices['C']}:")
    print("   Using fresh media is a critical part of quality control, but it is less likely to be the single solution than addressing the organism's growth rate.\n")
    
    # D
    print(f"D. {choices['D']}:")
    print(f"   This is a highly plausible solution. {facts['Suspected Pathogen']} {textwrap.fill(facts['Pathogen Characteristics'], 63)} The lab checked at {facts['Lab 1 Action']}. Standard procedure often calls for incubation for up to 72 hours. A longer incubation period would give the slow-growing Campylobacter a chance to form visible colonies.\n")

    # E
    print(f"E. {choices['E']}:")
    print("   This improves the statistical chance of isolation if the organism is rare, but it won't help if the colonies are simply too small to be seen at 48 hours. Longer incubation is more direct.\n")
    
    print("Conclusion:")
    print("-" * 35)
    print("The most direct and effective action to recover a potentially slow-growing pathogen from the existing culture plates would have been to allow more time for growth.")

# Execute the analysis
solve_microbiology_case()