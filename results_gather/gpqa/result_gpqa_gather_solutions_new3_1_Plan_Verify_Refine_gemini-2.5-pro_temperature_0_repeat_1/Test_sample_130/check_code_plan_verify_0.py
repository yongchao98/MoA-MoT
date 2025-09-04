def check_diels_alder_noesy():
    """
    This function checks the correctness of the answer to the Diels-Alder NOESY problem.
    It encapsulates the chemical logic required to solve the problem.
    
    1.  Identifies reactants and reaction type.
    2.  Assigns NMR signals to specific protons in the product.
    3.  Determines the major product (endo vs. exo) based on steric hindrance.
    4.  Analyzes spatial proximities in the major vs. minor product.
    5.  Identifies the unique NOESY cross-peak for the major product.
    6.  Compares the derived correct option with the provided answer.
    """
    
    # Step 1: Define problem parameters and signal assignments from the question
    signals = {
        "anhydride_H": {"desc": "a 2H singlet at ~3.5 ppm"},
        "vinylic_Me_H": {"desc": "a 6H singlet at ~1.7 ppm"},
        "bridgehead_Me_H": {"desc": "a 6H singlet at ~1.0 ppm"},
        "bridge_H": {"desc": "a 1H doublet at ~1.5 ppm"}
    }

    options = {
        "A": {"signal1": signals["bridgehead_Me_H"], "signal2": signals["bridge_H"]},
        "B": {"signal1": signals["bridgehead_Me_H"], "signal2": signals["vinylic_Me_H"]},
        "C": {"signal1": signals["bridge_H"], "signal2": signals["anhydride_H"]},
        "D": {"signal1": signals["vinylic_Me_H"], "signal2": signals["anhydride_H"]}
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer = "D"

    # Step 2: Analyze stereochemistry and determine the major product
    # The standard "Alder endo rule" predicts the 'endo' product is major due to electronic effects.
    # However, the diene (1,2,3,4-tetramethyl-1,3-cyclopentadiene) is extremely bulky.
    # The four methyl groups create severe steric hindrance for the 'endo' approach.
    # In such cases, steric hindrance can override the electronic preference, making the 'exo' product the major product.
    major_product_isomer = "exo"
    
    # Step 3: Analyze spatial proximities (NOESY interactions) for each isomer
    # In the 'exo' isomer, the anhydride ring is tucked under the C=C double bond.
    # This brings the anhydride protons very close to the vinylic methyl groups.
    exo_unique_interaction = {"anhydride_H", "vinylic_Me_H"}
    
    # In the 'endo' isomer, the anhydride ring is on the opposite side of the C=C double bond.
    # The anhydride protons are far from the vinylic methyls, but close to the C7 bridge protons.
    endo_unique_interaction = {"anhydride_H", "bridge_H"}

    # Step 4: Identify the expected cross-peak in the major product
    if major_product_isomer == "exo":
        expected_interaction = exo_unique_interaction
    else: # major_product_isomer == "endo"
        expected_interaction = endo_unique_interaction

    # Step 5: Determine the correct option based on the chemical analysis
    correct_option = None
    for option_key, option_signals in options.items():
        # Find the names of the proton groups corresponding to the signals in the option
        signal1_name = [k for k, v in signals.items() if v == option_signals["signal1"]][0]
        signal2_name = [k for k, v in signals.items() if v == option_signals["signal2"]][0]
        
        if {signal1_name, signal2_name} == expected_interaction:
            correct_option = option_key
            break

    # Step 6: Compare the LLM's answer with the derived correct answer
    if llm_answer == correct_option:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect. The correct answer is '{correct_option}'.\n"
        reason += "Reasoning:\n"
        reason += "1. The reaction is a Diels-Alder cycloaddition. Due to the severe steric hindrance of the four methyl groups on the diene, the 'exo' adduct is the major product, not the 'endo' adduct.\n"
        reason += "2. The question asks for a NOESY cross-peak present in the major product but absent in the minor. We must therefore analyze the structure of the major 'exo' product.\n"
        reason += "3. In the 'exo' adduct, the anhydride protons (2H singlet at ~3.5 ppm) are spatially close to the vinylic methyl groups (6H singlet at ~1.7 ppm).\n"
        reason += "4. In the minor 'endo' adduct, these two groups are far apart. Therefore, the unique cross-peak connects the signals for the anhydride protons and the vinylic methyl protons.\n"
        reason += f"5. This corresponds to option {correct_option}, which describes an interaction between '{signals[list(expected_interaction)[0]]['desc']}' and '{signals[list(expected_interaction)[1]]['desc']}'.\n"
        reason += f"The selected answer '{llm_answer}' describes an incorrect interaction."
        if llm_answer == "C":
            reason += " This interaction would be correct only if the 'endo' product were major, which is unlikely due to sterics."
        return reason

# Execute the check
result = check_diels_alder_noesy()
print(result)