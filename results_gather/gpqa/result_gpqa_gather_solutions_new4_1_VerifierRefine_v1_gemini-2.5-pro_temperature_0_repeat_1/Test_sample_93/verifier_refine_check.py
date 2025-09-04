def check_diels_alder_synthesis():
    """
    Checks the correctness of the provided answer for a chemical synthesis question.
    The script encodes the logic of the Diels-Alder reaction to verify the outcome.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_answer = 'C'

    # --- Step 1: Define the properties of the target molecule ---
    target_molecule = {
        "name": "methyl 2-propyl-1,2,4a,5,6,7,8,8a-octahydronaphthalene-1-carboxylate",
        "skeleton": "fused 6,6-ring system (decalin/octahydronaphthalene)",
        "double_bond_position": "C3-C4",
        "substituents": {
            "C1": "methyl carboxylate (-COOCH3)",
            "C2": "propyl (-CH2CH2CH3)"
        },
        "key_connectivity": "Substituents are on adjacent saturated carbons (C1, C2). The propyl group at C2 is adjacent to the C3=C4 double bond."
    }

    # --- Step 2: Analyze each option based on chemical principles ---

    # Analysis for Option A
    def analyze_A():
        analysis = {
            "starting_materials": "Cyclohexene and methyl 2,3-dimethylenehexanoate",
            "reaction_type": "Intermolecular Diels-Alder",
            "predicted_skeleton": "Spirocyclic compound (spiro[5.5]undecene derivative)",
            "reasoning": "The reaction of an exocyclic diene with a cyclic dienophile forms a spiro compound, not a fused ring system like the target.",
            "is_match": False
        }
        return analysis

    # Analysis for Option B
    def analyze_B():
        analysis = {
            "starting_materials": "1-vinylcyclohex-1-ene and methyl hex-2-ynoate",
            "reaction_type": "Intermolecular Diels-Alder",
            "predicted_product": "A tetrahydronaphthalene derivative (contains a cyclohexadiene ring with two double bonds).",
            "reasoning": "A Diels-Alder reaction with an alkyne dienophile produces a product with two double bonds in the newly formed ring. The target is an octahydronaphthalene with only one double bond.",
            "is_match": False
        }
        return analysis

    # Analysis for Option C
    def analyze_C():
        # Precursor: MeOOC-CH(2)=CH(3)-(CH2)4-CH(8)=CH(9)-CH(10)=CH(11)-Propyl
        analysis = {
            "starting_materials": "methyl (2E,8E,10E)-tetradeca-2,8,10-trienoate",
            "reaction_type": "Intramolecular Diels-Alder (IMDA)",
            "dienophile": "C2=C3 (activated by -COOCH3)",
            "diene": "C8=C9-C10=C11",
            "tether": "4-carbon chain (C4-C7), which is the correct length for a 6,6-fused system.",
            "cyclization": "New bonds form between C2-C11 and C3-C8.",
            "predicted_skeleton": "fused 6,6-ring system (decalin/octahydronaphthalene)",
            "predicted_double_bond": "Forms between C9 and C10 of the precursor.",
            "predicted_substituents": "'-COOCH3' from C2 and 'propyl' from C11 are on adjacent carbons.",
            "mapping_to_target": {
                "Product C1": "Precursor C2 (with -COOCH3)",
                "Product C2": "Precursor C11 (with propyl)",
                "Product C3=C4": "Precursor C10=C9"
            },
            "reasoning": "The mapping shows that the product has the -COOCH3 group at C1, the propyl group at C2, and the double bond at C3=C4. This perfectly matches all features of the target molecule.",
            "is_match": True
        }
        return analysis

    # Analysis for Option D
    def analyze_D():
        # Precursor: MeOOC-CH=CH-CH=CH-(CH2)4-CH=CH-Propyl
        analysis = {
            "starting_materials": "methyl (2E,4E,10Z)-tetradeca-2,4,10-trienoate",
            "reaction_type": "Intramolecular Diels-Alder (IMDA)",
            "diene": "C2=C3-C4=C5",
            "dienophile": "C10=C11",
            "cyclization": "New bonds form between C2-C11 and C5-C10.",
            "predicted_double_bond": "Forms between C3 and C4 of the precursor.",
            "reasoning": "The product's connectivity would be (fusion)-C(propyl)-C(COOMe)-C3=C4-(fusion). In the target molecule, the propyl group is adjacent to the double bond, not the ester group. This option produces a constitutional isomer, not the target molecule.",
            "is_match": False
        }
        return analysis

    # --- Step 3: Verify the LLM's answer ---
    analyses = {
        'A': analyze_A(),
        'B': analyze_B(),
        'C': analyze_C(),
        'D': analyze_D()
    }

    # Check if the provided answer is valid
    if llm_answer not in analyses:
        return f"Invalid answer option '{llm_answer}'. Please choose from A, B, C, or D."

    chosen_analysis = analyses[llm_answer]

    if chosen_analysis["is_match"]:
        return "Correct"
    else:
        # Find the correct answer to provide a more helpful message
        correct_option = None
        for option, analysis in analyses.items():
            if analysis["is_match"]:
                correct_option = option
                break
        
        error_reason = f"The provided answer '{llm_answer}' is incorrect. "
        error_reason += f"Reasoning: {chosen_analysis['reasoning']}. "
        
        if correct_option:
            correct_analysis = analyses[correct_option]
            error_reason += f"The correct answer is '{correct_option}' because its intramolecular Diels-Alder reaction directly yields the target molecule with the correct skeleton, substituent positions, and double bond location, as explained: {correct_analysis['reasoning']}."
        else:
            error_reason += "Based on standard reaction analysis, none of the options appear to correctly produce the target molecule."
            
        return error_reason

# Run the check and print the result
result = check_diels_alder_synthesis()
print(result)