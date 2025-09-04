import re

def check_organic_reaction_answer(question, llm_answer):
    """
    Checks the correctness of an organic chemistry reaction sequence explanation.

    This function does not simulate the reaction but validates the chemical
    reasoning presented in the explanation.
    """
    
    # Extract the explanation and the final choice (e.g., <<<C>>>)
    try:
        explanation = llm_answer.split("<<<")[0]
        final_choice = llm_answer.split("<<<")[1].split(">>>")[0].strip()
    except IndexError:
        return "Incorrect. The answer format is invalid. It should end with '<<<X>>>'."

    errors = []

    # --- Check 1: Reaction Step Plausibility ---
    # Step 1: H2O -> A (SN1 Solvolysis)
    step1_keywords = ["solvolysis", "sn1", "alcohol"]
    if not all(re.search(kw, explanation, re.IGNORECASE) for kw in step1_keywords):
        errors.append("Step 1 explanation is incomplete. It should describe an SN1/solvolysis reaction forming an alcohol.")

    # Step 2: A + PDC -> B (Oxidation)
    step2_keywords = ["pdc", "oxidizing agent", "oxidizes", "secondary alcohol", "ketone"]
    if not all(re.search(kw, explanation, re.IGNORECASE) for kw in step2_keywords):
        errors.append("Step 2 explanation is incomplete. It should describe the oxidation of a secondary alcohol to a ketone using PDC.")

    # Step 3: B + H2CPPh3 -> C (Wittig Reaction)
    step3_keywords = ["wittig", "ketone", "alkene", "h₂cpph₃|methylenetriphenylphosphorane"]
    if not all(re.search(kw, explanation, re.IGNORECASE) for kw in step3_keywords):
        errors.append("Step 3 explanation is incomplete. It should describe a Wittig reaction converting a ketone to an alkene.")

    # --- Check 2: Rearrangement Mechanism Logic ---
    # Step 4: C + TsOH -> D (Acid-catalyzed rearrangement)
    step4_keywords = ["tsoh|p-toluenesulfonic acid", "acid", "protonates", "carbocation", "rearrangement", "methyl shift", "wagner-meerwein"]
    if not all(re.search(kw, explanation, re.IGNORECASE) for kw in step4_keywords):
        errors.append("Step 4 explanation is incomplete. It should detail an acid-catalyzed rearrangement involving a carbocation and a Wagner-Meerwein (methyl) shift.")

    # --- Check 3: Substituent Transformation Consistency ---
    # The core logic is the transformation of substituents.
    # Start: 3a,4a-dimethyl
    # End: 3a,5,5-trimethyl
    start_pattern = "3a,4a-dimethyl"
    end_pattern = "3a,5,5-trimethyl"
    
    if not (re.search(start_pattern, explanation, re.IGNORECASE) and re.search(end_pattern, explanation, re.IGNORECASE)):
        errors.append(f"The explanation fails to mention the starting ({start_pattern}) and final ({end_pattern}) substituent patterns.")
    
    # Check if the logic for the change is explained
    shift_logic = "methyl group from c4a migrates to c5|methyl group from c4a has moved to c5"
    wittig_methyl_logic = "former =ch₂ group has become a new methyl group"
    if not (re.search(shift_logic, explanation, re.IGNORECASE) and re.search(wittig_methyl_logic, explanation, re.IGNORECASE)):
        errors.append("The explanation does not correctly justify the formation of the 5,5-dimethyl group via methyl shift and protonation of the methylene group.")

    # --- Check 4: Final Answer Coherence ---
    # The explanation concludes that the product is named '3a,5,5-trimethyl...' and corresponds to option C.
    final_name_in_text = "3a,5,5-trimethyl-1,2,3,3a,5,6,7,8-octahydrocyclopenta[1,4]cyclobuta[1,2]benzene"
    if final_name_in_text.lower() not in explanation.lower():
        errors.append("The explanation does not conclude with the specific name of the final product.")
    
    if final_choice != 'C':
        errors.append(f"The final choice is '{final_choice}', but the explanation points to the structure corresponding to option C.")

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. Reason(s):\n- " + "\n- ".join(errors)


# The user's question and the provided LLM answer
question = """
Identify product D in the following series of reactions.
5-bromo-3a,4a-dimethyldecahydrocyclopenta[1,4]cyclobuta[1,2]benzene + H2O ---> A
A + PDC ---> B
B + H2CPPh3 ---> C
C + TsOH ---> D

A) 3a,4,5a-trimethyl-1,2,3,3a,5a,6,7,8-octahydrocyclopenta[c]pentalene
B) 3a,4a,5,5-tetramethyl-2,3,3a,4,4a,5-hexahydro-1H-cyclobuta[1,2:1,4]di[5]annulene
C) 3a,5,5-trimethyl-1,2,3,3a,5,6,7,8-octahydrocyclopenta[1,4]cyclobuta[1,2]benzene
D) 3a,5-dimethyldecahydrocyclopenta[1,4]cyclobuta[1,2]benzene
"""

llm_answer = """
An explanation of the reaction sequence:

1.  **Starting Material to A:** The starting material is a tricyclic alkyl halide. Reacting it with water (H₂O), a weak nucleophile and polar protic solvent, leads to a solvolysis reaction, specifically an SN1 reaction. The bromide (Br⁻) acts as a good leaving group, and it is replaced by a hydroxyl group (-OH) to form a secondary alcohol, **A**.

2.  **A to B:** Pyridinium dichromate (PDC) is a mild oxidizing agent. It selectively oxidizes the secondary alcohol (**A**) to a ketone (**B**) without affecting other parts of the molecule.

3.  **B to C:** The reaction of the ketone (**B**) with methylenetriphenylphosphorane (H₂C=PPh₃) is a Wittig reaction. This reaction converts the carbonyl group (C=O) into an alkene. Specifically, the oxygen atom is replaced by the =CH₂ group, forming an exo-methylene group (a terminal alkene), which is product **C**.

4.  **C to D:** The final step involves treating the alkene (**C**) with a strong acid, p-toluenesulfonic acid (TsOH). This initiates a carbocation-mediated rearrangement:
    *   **Protonation:** The acid protonates the terminal alkene at the less substituted carbon (the CH₂ group), forming the more stable tertiary carbocation at C5. The former =CH₂ group is now a methyl group.
    *   **Wagner-Meerwein Rearrangement:** The molecule now has a tertiary carbocation at C5, which is adjacent to a quaternary carbon at C4a (bearing a methyl group). To form an even more stable tertiary carbocation and relieve steric strain, a 1,2-methyl shift occurs. The methyl group from C4a migrates to C5.
    *   **Result of Shift:** This migration creates a new, stable tertiary carbocation at C4a and results in a gem-dimethyl group (two methyl groups on the same carbon) at C5.
    *   **Elimination:** Finally, a proton is eliminated from a carbon adjacent to the C4a carbocation to form a stable, substituted internal alkene.

The key structural change is the substituent pattern. The starting material was `3a,4a-dimethyl`. After the rearrangement, the methyl group from C4a has moved to C5, and the former =CH₂ group has become a new methyl group at C5. The final substituent pattern is `3a-methyl` and `5,5-dimethyl`. This corresponds to the name **3a,5,5-trimethyl-1,2,3,3a,5,6,7,8-octahydrocyclopenta[1,4]cyclobuta[1,2]benzene**.

Therefore, product **D** is option C.
<<<C>>>
"""

# Run the check
result = check_organic_reaction_answer(question, llm_answer)
print(result)