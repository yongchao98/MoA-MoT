import re

def check_analytical_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the answer to the analytical chemistry question.

    The function models the properties of each option (precursor, enantiomers, etc.)
    and verifies if the chosen answer's properties match the experimental observations
    from NMR and LC-MS.
    """

    # Step 1: Define the properties of each possible explanation based on chemical principles.
    # The options are taken directly from the question text.
    options_properties = {
        'A': {
            "name": "Contamination with a precursor",
            "is_isomer": False,  # Precursors have a different (usually lower) mass.
            "is_separable_by_achiral_lc": True,
            "is_distinguishable_by_achiral_nmr": True
        },
        'B': {
            "name": "Mixture of enantiomers",
            "is_isomer": True,   # Enantiomers are isomers (same mass).
            "is_separable_by_achiral_lc": False, # Not separable by standard (achiral) LC.
            "is_distinguishable_by_achiral_nmr": False # Not distinguishable by standard (achiral) NMR.
        },
        'C': {
            "name": "'Double coupling' side product",
            "is_isomer": False,  # This product has a different (higher) mass.
            "is_separable_by_achiral_lc": True,
            "is_distinguishable_by_achiral_nmr": True
        },
        'D': {
            "name": "Mixture of diastereoisomers",
            "is_isomer": True,   # Diastereoisomers are isomers (same mass).
            "is_separable_by_achiral_lc": True,  # Separable by standard (achiral) LC.
            "is_distinguishable_by_achiral_nmr": True # Distinguishable by standard (achiral) NMR.
        }
    }

    # Step 2: Define the constraints derived from the experimental observations in the question.
    experimental_constraints = {
        "is_isomer": {
            "value": True,
            "reason": "LC-MS analysis shows both peaks have the same mass spectrum, consistent with the expected molecule."
        },
        "is_separable_by_achiral_lc": {
            "value": True,
            "reason": "LC-MS analysis shows two clearly defined peaks."
        },
        "is_distinguishable_by_achiral_nmr": {
            "value": True,
            "reason": "1H NMR spectrum shows two peaks for the same alpha-proton (and spin-spin coupling is ruled out)."
        }
    }

    # Step 3: Extract the final answer from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not find the final answer in the required format, e.g., <<<D>>>."
    
    final_answer_key = match.group(1)
    
    if final_answer_key not in options_properties:
        return f"Failure: The extracted answer '{final_answer_key}' is not a valid option (A, B, C, or D)."

    chosen_option_props = options_properties[final_answer_key]
    chosen_option_name = chosen_option_props['name']

    # Step 4: Check if the chosen answer satisfies all experimental constraints.
    for constraint_key, constraint_data in experimental_constraints.items():
        expected_value = constraint_data["value"]
        reason_for_constraint = constraint_data["reason"]
        
        if chosen_option_props[constraint_key] != expected_value:
            return (f"Incorrect. The chosen answer '{final_answer_key}: {chosen_option_name}' is wrong.\n"
                    f"Constraint not satisfied: The experimental data shows that '{reason_for_constraint}', which means the correct explanation must have the property '{constraint_key}' as {expected_value}.\n"
                    f"Reason for failure: A '{chosen_option_name}' has the property '{constraint_key}' as {chosen_option_props[constraint_key]}, which contradicts the observation.")

    # Step 5: Verify that only one option satisfies all constraints to ensure the question is well-posed.
    correct_options = []
    for key, props in options_properties.items():
        is_match = all(props[k] == experimental_constraints[k]['value'] for k in experimental_constraints)
        if is_match:
            correct_options.append(key)
    
    if len(correct_options) > 1:
        return f"Ambiguous Question: The chosen answer '{final_answer_key}' fits all criteria, but so does/do option(s) {', '.join(o for o in correct_options if o != final_answer_key)}. The question may be flawed."
    
    if final_answer_key not in correct_options:
        # This case should be caught by the loop in step 4, but it's a good final check.
        return f"Incorrect. The chosen answer '{final_answer_key}' does not fit the experimental criteria. The only correct answer is {correct_options[0]}."

    # If all checks pass, the answer is correct.
    return "Correct"

# The provided answer to be checked
llm_answer = """
Here is a step-by-step analysis of the problem, incorporating careful points derived from the experimental evidence.

### Step 1: Analyze the Mass Spectrometry (MS) Data

This is the most definitive piece of evidence and should be the first point of analysis.
*   **Observation:** The LC-MS analysis shows two distinct peaks in the chromatogram, but both peaks yield the *same mass spectrum*. This mass is also consistent with the expected molecule.
*   **Careful Point:** This observation proves that the two chemical species are **isomers**. They have the same chemical formula and molecular weight but a different three-dimensional arrangement of atoms. This single fact allows us to immediately eliminate any options that involve a change in mass.
    *   **Eliminates A) The compound is contaminated with a precursor.** A precursor is an intermediate from the synthesis and would have a different (usually lower) molecular weight.
    *   **Eliminates C) 'Double coupling' has occurred.** This side reaction would result in a product with an extra amino acid, leading to a different (higher) molecular weight.

### Step 2: Analyze the Chromatography (LC) and NMR Data

With precursors and mass-variant side products ruled out, we must distinguish between the two remaining types of isomers: enantiomers (B) and diastereoisomers (D).
*   **Observation:** The two species are separated by liquid chromatography (LC), appearing as two distinct peaks. The ¹H NMR spectrum shows two distinct peaks for a single alpha-proton.
*   **Careful Point:** These observations prove that the two isomers are physically and chemically distinct under standard (achiral) analytical conditions.
    *   **Evaluate B) The crude compound exists as a mixture of enantiomers:** Enantiomers are non-superimposable mirror images. In an achiral environment (like a standard LC column or a standard NMR solvent), they have identical physical and chemical properties. Therefore, they would **not** be separated by a standard LC column and would **not** be distinguishable in a standard NMR experiment. The experimental data directly contradicts this.
    *   **Evaluate D) The crude compound exists as a mixture of diastereoisomers:** Diastereoisomers are stereoisomers that are *not* mirror images of each other. They have different physical properties and are chemically distinct. This means they can be separated by standard LC, and their corresponding protons will exist in different chemical environments, leading to distinct signals in an NMR spectrum. This explanation is perfectly consistent with all the provided data.

### Step 3: Synthesize the Conclusion

By systematically applying the constraints from the experimental data, we can confidently arrive at the correct answer.
1.  The identical mass of the two species confirms they are isomers, ruling out precursors and side-products.
2.  The fact that they are separable by LC and distinguishable by NMR proves they are not enantiomers.
3.  The only remaining possibility that fits all the evidence is a mixture of diastereoisomers. This is a common phenomenon in peptide synthesis, often arising from either epimerization (inversion of a single chiral center) or the presence of stable *cis/trans* isomers (rotamers) around a peptide bond.

<<<D>>>
"""

# Run the check
result = check_analytical_chemistry_answer(llm_answer)
print(result)