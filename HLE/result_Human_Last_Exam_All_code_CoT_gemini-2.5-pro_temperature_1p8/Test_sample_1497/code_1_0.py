def solve_violin_puzzle():
    """
    Solves the violin timbre puzzle by mapping physical variations
    to the given parameter groups.
    """

    # Define the parameter groups and their primary physical meaning.
    parameter_groups = {
        'i': "F (fundamental frequency of the open string)",
        'ii': "ν (controls strength of the string's higher harmonics)",
        'iii': "a_1, f_1 (parameters for the first main body resonance)",
        'iv': "μ, a_2, f_2 (parameters for the second/higher body resonances)"
    }

    # Analyze each variation and determine the corresponding parameter group.
    
    # (1) Sul ponticello enhances string harmonics -> ν -> group ii
    map_1 = 'ii'

    # (2) Mute dampens/mass-loads body resonances -> a_1, f_1 -> group iii
    map_2 = 'iii'
    # Adding mass to the bridge lowers its resonant frequency f_1.
    direction_2 = 'down'

    # (3) Helium changes the speed of sound, altering air resonances -> f_m -> group iv
    # This affects higher-frequency air modes significantly.
    map_3 = 'iv'

    # (4) Changing from A-string to E-string changes the open string fundamental -> F -> group i
    map_4 = 'i'

    # Assemble the final answer in the specified format.
    final_answer = f"{map_1},{map_2},{map_3},{map_4},{direction_2}"

    # Print the reasoning and the final answer.
    print("This script analyzes how four variations in violin playing affect its timbre,")
    print("mapping each variation to a group of parameters from the provided waveform model.")
    print("-" * 70)
    print(f"1. 'sul ponticello': Primarily affects {parameter_groups[map_1]}. Maps to group '{map_1}'.")
    print(f"2. Bridge Mute: Primarily affects {parameter_groups[map_2]}. Maps to group '{map_2}'.")
    print(f"   The direction of change for f_1 (resonant frequency) is '{direction_2}'.")
    print(f"3. Helium Environment: Primarily affects {parameter_groups[map_3]}. Maps to group '{map_3}'.")
    print(f"4. Playing on E-string: Primarily affects {parameter_groups[map_4]}. Maps to group '{map_4}'.")
    print("-" * 70)
    print("The final answer in the specified format is:")
    print(final_answer)

solve_violin_puzzle()