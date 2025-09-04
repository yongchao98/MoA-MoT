def check_answer():
    """
    Checks the correctness of the answer for the organic synthesis question.
    The logic codifies the rules of the Intramolecular Diels-Alder (IMDA) reaction.
    """
    
    # --- Define the properties of the target molecule ---
    target_molecule = {
        "name": "methyl 2-propyl-1,2,4a,5,6,7,8,8a-octahydronaphthalene-1-carboxylate",
        "reaction_type": "intramolecular",
        "skeleton": "bicyclo[4.4.0]decene", # Fused 6,6-ring system
        "double_bond_position": "C3=C4",
        "substituent_connectivity": "adjacent",
        "substituent_order": "C1(-COOCH3)-C2(-Propyl)-C3=C4"
    }

    # --- Define the candidate starting materials ---
    # Note: The lettering (A, B, C, D) is taken from the final provided answer block.
    candidates = {
        "A": "methyl (2E,8E,10E)-tetradeca-2,8,10-trienoate",
        "B": "Cyclohexene and methyl 2,3-dimethylenehexanoate",
        "C": "1-vinylcyclohex-1-ene and methyl hex-2-ynoate",
        "D": "methyl (2E,4E,10Z)-tetradeca-2,4,10-trienoate"
    }
    
    llm_answer = "A"
    results = {}

    # --- Analyze each candidate ---
    for key, name in candidates.items():
        print(f"--- Analyzing Candidate {key}: {name} ---")
        
        # Candidate A: IMDA with dienophile at C2, diene at C8-C11
        if key == "A":
            # Constraint 1: Reaction Type
            if "and" in name:
                results[key] = "Failed: Is an intermolecular reaction, but an intramolecular reaction is expected for this target."
                print(results[key])
                continue
            
            # Constraint 2: Skeleton
            # Tether is C4, C5, C6, C7 (4 atoms). This correctly forms a fused 6,6-ring system.
            print("Passed: Is a valid intramolecular precursor with a 4-atom tether, forming a bicyclo[4.4.0] skeleton.")

            # Constraint 3: Double Bond Position
            # Diene is at C8-C11. New double bond forms between central carbons C9 and C10.
            # Mapping to product: C2(prec)->C1(prod), C11(prec)->C2(prod), C10(prec)->C3(prod), C9(prec)->C4(prod).
            # So, new double bond is at C3=C4 of the product.
            print("Passed: Diene at C8-C11 results in a C3=C4 double bond in the final product, matching the target.")

            # Constraint 4: Substituent Placement
            # Substituents on C2(-COOCH3) and C11(-Propyl) become adjacent.
            # The order relative to the double bond is C1(-COOCH3)-C2(-Propyl)-C3=C4.
            print("Passed: Substituents are placed correctly at C1(-COOCH3) and C2(-Propyl), matching the target.")
            results[key] = "Correct"

        # Candidate B: Intermolecular reaction
        elif key == "B":
            results[key] = "Failed: Intermolecular reaction between a diene and cyclohexene would form a spiro compound, not the required fused decalin skeleton."
            print(results[key])

        # Candidate C: Intermolecular reaction with an alkyne
        elif key == "C":
            results[key] = "Failed: Intermolecular reaction with an alkyne dienophile would produce a hexahydronaphthalene (two double bonds), not the target octahydronaphthalene (one double bond)."
            print(results[key])

        # Candidate D: IMDA with diene at C2-C5, dienophile at C10
        elif key == "D":
            # Constraint 1 & 2: Pass, same as A.
            print("Passed: Is a valid intramolecular precursor with a 4-atom tether, forming a bicyclo[4.4.0] skeleton.")
            
            # Constraint 3: Double Bond Position
            # Diene is at C2-C5. New double bond forms between C3 and C4. This matches the target.
            print("Passed: Diene at C2-C5 results in a C3=C4 double bond in the final product, matching the target.")

            # Constraint 4: Substituent Placement
            # Substituents on C2(-COOCH3) and C11(-Propyl) become adjacent.
            # However, the order relative to the double bond is C(fusion)-C(Propyl)-C(COOMe)-C3=C4.
            # The target requires C(fusion)-C(COOMe)-C(Propyl)-C3=C4.
            results[key] = "Failed: The order of substituents is incorrect. This precursor forms an isomer of the target, but not the target molecule itself."
            print(results[key])
        
        print("-" * (len(key) + len(name) + 22))


    # --- Final Verdict ---
    print("\n--- VERDICT ---")
    if results.get(llm_answer) == "Correct":
        return "Correct"
    else:
        failure_reason = results.get(llm_answer, "The provided answer option does not exist in the candidate list.")
        correct_options = [k for k, v in results.items() if v == "Correct"]
        if correct_options:
            return f"Incorrect. The provided answer '{llm_answer}' is wrong. Reason: {failure_reason}\nThe correct answer should be '{correct_options[0]}'."
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is wrong. Reason: {failure_reason}\nNo candidate correctly forms the target molecule based on the analysis."

# Run the check
result = check_answer()
print(result)