import textwrap

def analyze_lab_scenario():
    """
    Analyzes a microbiology case study to identify the likely process failure
    and determine the best corrective action from a list of options.
    """
    # The problem describes a scenario where a lab tries to isolate Campylobacter.
    # Key details from the prompt:
    # 1. Symptom: Bloody diarrhea (suggests Campylobacter).
    # 2. Media: Campy-Cefex agar (selective for Campylobacter).
    # 3. Incubation: Microaerophilic container at 42 degrees for 2 days (specific for Campylobacter).
    # 4. Result: Overgrowth of a contaminant (Bacillus) and failure to identify the pathogen.
    # 5. The core issue is that the selective medium failed to inhibit the contaminant.

    print("Step 1: Analyze the laboratory protocol and results.")
    print("The protocol used by the first lab was correct for isolating Campylobacter.")
    print(f"- It used selective media (Campy-Cefex agar).")
    print(f"- It used the correct incubation conditions, including the specific temperature of 42 degrees Celsius and a microaerophilic atmosphere.")
    print(f"- The incubation time of 2 days is standard.")
    print("However, the result was the growth of 'Bacillus species', a common contaminant, which overgrew the plate. This indicates the selective media did not perform its function of suppressing non-target bacteria.\n")

    print("Step 2: Evaluate each possible solution.\n")

    # Dictionary of choices and their analysis
    analysis = {
        'A. Obtain a fresh sample': 'This does not solve the problem of how to correctly process the original sample. It circumvents the issue rather than fixing the flawed procedure.',
        
        'B. Decrease sample processing time': 'While good practice, the primary problem identified was contaminant overgrowth, not the death of the target organism due to delay. This would have been a less impactful change.',
        
        'C. Used Fresh Campylobacter plates': 'This is the most probable solution. The antibiotics in selective media like Campy-Cefex agar can degrade over time if plates are old or improperly stored. Degraded antibiotics would fail to inhibit contaminant bacteria like Bacillus, allowing them to overgrow the slower-growing Campylobacter. Using fresh plates ensures the selective agents are at full potency.',
        
        'D. Incubated the sample for longer': 'This would worsen the problem. The fast-growing Bacillus contaminant already produced large colonies in 2 days. More time would only allow for more overgrowth, further obscuring any potential Campylobacter colonies.',
        
        'E. Increase the number of plates': 'If the entire batch of media is ineffective, plating on more plates from the same bad batch would not help. It would simply yield more overgrown plates.'
    }

    # Print the analysis for each option
    for option, explanation in analysis.items():
        print(f"--- Analysis of Option {option} ---")
        # textwrap helps format the explanation nicely
        wrapped_explanation = textwrap.fill(explanation, width=70)
        print(wrapped_explanation)
        print("-" * 35 + "\n")

    print("Step 3: Conclude the most likely corrective action.")
    print("The most direct cause of the isolation failure was the ineffectiveness of the selective agar, allowing a contaminant to flourish. Therefore, ensuring the potency of the media by using fresh plates would be the single most effective way the first laboratory could have recovered the organism.")

# Execute the analysis
analyze_lab_scenario()
<<<C>>>