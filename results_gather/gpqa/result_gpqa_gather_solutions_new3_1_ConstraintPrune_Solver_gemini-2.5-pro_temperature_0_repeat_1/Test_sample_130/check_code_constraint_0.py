def check_answer():
    """
    Checks the correctness of the final answer for the given chemistry problem.

    The logic is based on a step-by-step chemical analysis:
    1.  Identify reactants and products.
    2.  Assign proton NMR signals.
    3.  Determine the major product (endo vs. exo). This is the key step.
    4.  Identify the unique NOESY cross-peak in the major product.
    5.  Compare the expected cross-peak with the chosen answer.
    """
    
    # The final answer provided by the analysis
    final_answer_letter = "C"
    
    # --- Chemical Knowledge Base ---
    
    # 1. NMR Signal Assignments
    signal_assignments = {
        "H_anhydride": "2H singlet at ~3.5 ppm",
        "Me_vinylic": "6H singlet at ~1.7 ppm",
        "Me_bridgehead": "6H singlet at ~1.0 ppm",
        "H_bridge": "1H doublet at ~1.5 ppm"
    }

    # 2. Spatial Proximities (NOE) in each isomer
    # These are the unique interactions that distinguish the isomers.
    noe_in_exo_isomer = ("H_anhydride", "Me_vinylic")
    noe_in_endo_isomer = ("H_anhydride", "H_bridge") # Also with Me_bridgehead, but this is a key one.

    # 3. Major Product Determination
    # Due to severe steric hindrance from the four methyl groups on the diene,
    # the exo transition state is favored over the endo.
    major_product_isomer = "exo"

    # --- Verification Logic ---

    # Step 1: Determine the expected unique NOE based on the major product
    if major_product_isomer == "exo":
        expected_interaction = noe_in_exo_isomer
        reasoning_major_product = "The major product is the EXO adduct due to severe steric hindrance from the tetramethyl-substituted diene."
        reasoning_noe = f"In the EXO adduct, the anhydride protons ({signal_assignments['H_anhydride']}) are spatially close to the vinylic methyl groups ({signal_assignments['Me_vinylic']})."
    else: # (Hypothetical case if endo were major)
        expected_interaction = noe_in_endo_isomer
        reasoning_major_product = "The major product is the ENDO adduct (standard Alder rule)."
        reasoning_noe = f"In the ENDO adduct, the anhydride protons ({signal_assignments['H_anhydride']}) are spatially close to the bridge protons ({signal_assignments['H_bridge']})."

    # Step 2: Get the signal descriptions for the expected interaction
    proton1_signal = signal_assignments[expected_interaction[0]]
    proton2_signal = signal_assignments[expected_interaction[1]]

    # Step 3: Define the options from the question
    options = {
        "A": "A 6H singlet at ~1 ppm and a 6H singlet at ~1.7 ppm",
        "B": "A 1H doublet at ~1.5 ppm and a 2H singlet at ~3.5 ppm",
        "C": "A 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm",
        "D": "A 6H singlet at ~1 ppm and a 1H doublet at ~1.5 ppm"
    }
    
    chosen_answer_text = options.get(final_answer_letter)

    # Step 4: Check if the chosen answer matches the expected interaction
    # We check if both expected signal descriptions are present in the chosen answer text.
    is_correct = (proton1_signal.split(" at ")[0] in chosen_answer_text and
                  proton2_signal.split(" at ")[0] in chosen_answer_text)

    if is_correct:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The provided answer '{final_answer_letter}' is wrong.\n"
            f"Reasoning:\n"
            f"1. {reasoning_major_product}\n"
            f"2. {reasoning_noe}\n"
            f"3. Therefore, the unique cross-peak in the major product should connect the signals for '{proton1_signal}' and '{proton2_signal}'.\n"
            f"4. The chosen answer '{chosen_answer_text}' does not describe this interaction."
        )
        return error_message

# Run the check
result = check_answer()
print(result)