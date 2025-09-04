import re

def check_chemistry_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the given organic chemistry problem.

    The function verifies the answer based on two key chemical principles:
    1.  **Isomerism in Reaction 1:** The product of a rearrangement reaction (oxy-Cope) must be an isomer of the starting material. This is checked by comparing molecular formulas.
    2.  **Product Form in Reaction 2:** The product of the Ireland-Claisen rearrangement under basic conditions (LDA) without an acidic workup must be the corresponding salt, not the free acid.

    Args:
        final_answer_text: The text of the final answer, which includes the reasoning and the final choice in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a detailed reason for why the answer is incorrect.
    """
    # --- 1. Define Chemical Data and Constraints ---

    # Molecular formulas are determined from chemical structure analysis.
    # This is a key constraint for Reaction 1 (rearrangement must conserve atoms).
    molecular_formulas = {
        "1-vinylspiro[3.5]non-5-en-1-ol": "C11H16O",
        "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one": "C11H16O",
        "decahydro-7H-benzo[7]annulen-7-one": "C11H18O",
    }

    # Define the multiple-choice options as presented in the problem.
    options = {
        "A": {"A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one", "B": "lithium 3-ethylpent-4-enoate"},
        "B": {"A": "decahydro-7H-benzo[7]annulen-7-one", "B": "3-ethylpent-4-enoic acid"},
        "C": {"A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one", "B": "3-ethylpent-4-enoic acid"},
        "D": {"A": "decahydro-7H-benzo[7]annulen-7-one", "B": "lithium 3-ethylpent-4-enoate"},
    }

    # --- 2. Parse the Provided Final Answer ---

    # Extract the chosen option letter (e.g., 'A', 'B', 'C', or 'D').
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Failure: Could not parse the final answer choice from the text. Expected format is <<<X>>>."
    
    llm_choice = match.group(1)

    # --- 3. Determine the Correct Answer Programmatically ---

    # Constraint for Reaction 1: Product A must be an isomer of the starting material.
    start_material_formula = molecular_formulas["1-vinylspiro[3.5]non-5-en-1-ol"]
    correct_product_A = None
    for product_name, product_formula in molecular_formulas.items():
        if product_name != "1-vinylspiro[3.5]non-5-en-1-ol" and product_formula == start_material_formula:
            correct_product_A = product_name
            break
    
    if not correct_product_A:
        return "Checker Error: Could not determine the correct product A based on isomerism."

    # Constraint for Reaction 2: Product B must be a lithium salt due to the use of LDA and no acidic workup.
    # The product should be a carboxylate salt, not a carboxylic acid.
    correct_product_B = "lithium 3-ethylpent-4-enoate"

    # Find which option letter corresponds to the correct products.
    correct_option_letter = None
    for option, products in options.items():
        if products["A"] == correct_product_A and products["B"] == correct_product_B:
            correct_option_letter = option
            break

    if not correct_option_letter:
        return "Checker Error: No option matches the derived correct products."

    # --- 4. Compare and Return the Result ---

    if llm_choice == correct_option_letter:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer chose option {llm_choice}, but the correct option is {correct_option_letter}.\n\n"
            f"**Reasoning:**\n"
            f"1.  **Constraint for Product A:** The reaction is an isomerization (anionic oxy-Cope rearrangement), so the product must have the same molecular formula as the starting material (C11H16O). "
            f"'{correct_product_A}' has the formula C11H16O, while 'decahydro-7H-benzo[7]annulen-7-one' has the formula C11H18O. Therefore, Product A must be '{correct_product_A}'. This eliminates options B and D.\n"
            f"2.  **Constraint for Product B:** The reaction (Ireland-Claisen rearrangement) uses a strong lithium base (LDA) and does not specify an acidic workup. The product is therefore the deprotonated carboxylate salt. "
            f"Product B must be '{correct_product_B}', not the free acid '3-ethylpent-4-enoic acid'. This eliminates option C.\n"
            f"3.  **Conclusion:** The only option that satisfies both constraints is {correct_option_letter}."
        )
        return reason

# The final answer provided by the user to be checked.
final_answer_to_check = """
Here is a step-by-step analysis to determine the correct answer.

### **Analysis of Reaction 1: Synthesis of Product A**

1.  **Starting Material and Reagents:**
    *   **Reactant:** `1-vinylspiro[3.5]non-5-en-1-ol`. This molecule is a tertiary alcohol and, more importantly, a 3-hydroxy-1,5-diene system, which is the characteristic substrate for an oxy-Cope rearrangement.
    *   **Reagents:** `KH` (Potassium Hydride) is a strong base that deprotonates the alcohol to form a potassium alkoxide. `H+` indicates an acidic workup step.

2.  **Mechanism:**
    *   The formation of an alkoxide from a 3-hydroxy-1,5-diene under thermal or base-catalyzed conditions initiates an **anionic oxy-Cope rearrangement**. This is a type of [3,3]-sigmatropic rearrangement.
    *   The rearrangement is highly favorable as it breaks a bond within the strained four-membered cyclobutane ring.
    *   The initial rearrangement product is a 10-membered ring containing an enolate.
    *   Under the reaction conditions (especially with the acidic workup), this flexible 10-membered ring intermediate can undergo a subsequent **transannular cyclization** to form a more thermodynamically stable, bridged bicyclic system.
    *   The acidic workup (H+) protonates the final enolate, which tautomerizes to the stable ketone product.

3.  **Identifying Product A:**
    *   Let's check the molecular formulas. The starting material is **C₁₁H₁₆O**. Any rearrangement product must be an isomer and have the same formula.
    *   Option 1 for A: `(E)-bicyclo[5.3.1]undec-1(11)-en-4-one`. This bicyclic enone has the formula **C₁₁H₁₆O**. This is a valid isomer and is the known product of this tandem reaction sequence.
    *   Option 2 for A: `decahydro-7H-benzo[7]annulen-7-one`. This saturated bicyclic ketone has the formula **C₁₁H₁₈O**. This is not an isomer of the starting material (it has two extra hydrogens), so it cannot be the product.
    *   **Conclusion for A:** Product A must be **(E)-bicyclo[5.3.1]undec-1(11)-en-4-one**. This eliminates options B and D.

### **Analysis of Reaction 2: Synthesis of Product B**

1.  **Starting Material and Reagents:**
    *   **Reactant:** `(E)-pent-2-en-1-ol`, an allylic alcohol.
    *   **Reagents:** `acetyl bromide` (an acylating agent) and `LDA` (Lithium Diisopropylamide), a strong, non-nucleophilic base.

2.  **Mechanism:**
    *   This set of reagents is characteristic of the **Ireland-Claisen rearrangement**.
    *   The mechanism involves the *in situ* formation of an allylic ester, `(E)-pent-2-en-1-yl acetate`.
    *   The strong base, LDA, then deprotonates the α-carbon of the ester's acetyl group to form a lithium enolate.
    *   This lithium enolate undergoes a [3,3]-sigmatropic rearrangement.

3.  **Identifying Product B:**
    *   The direct product of the Ireland-Claisen rearrangement is a carboxylate salt.
    *   The problem specifies the base is LDA (a lithium base) and, crucially, **does not mention a final acidic workup step** for this second reaction.
    *   Therefore, the carboxylic acid product will remain in its deprotonated form as a lithium salt.
    *   The rearranged carbon skeleton corresponds to 3-ethylpent-4-enoic acid. As a salt, its name is 3-ethylpent-4-enoate.
    *   **Conclusion for B:** The product is the lithium salt, **lithium 3-ethylpent-4-enoate**. This rules out the options that list the free carboxylic acid.

### **Final Conclusion**

*   Product A is **(E)-bicyclo[5.3.1]undec-1(11)-en-4-one**.
*   Product B is **lithium 3-ethylpent-4-enoate**.

Matching these two products with the given options:

*   **A) A = (E)-bicyclo[5.3.1]undec-1(11)-en-4-one, B = lithium 3-ethylpent-4-enoate**
*   B) A = decahydro-7H-benzo[7]annulen-7-one, B = 3-ethylpent-4-enoic acid
*   C) A = (E)-bicyclo[5.3.1]undec-1(11)-en-4-one, B = 3-ethylpent-4-enoic acid
*   D) A = decahydro-7H-benzo[7]annulen-7-one, B = lithium 3-ethylpent-4-enoate

Only option A correctly identifies both products.

<<<A>>>
"""

# Run the check
result = check_chemistry_answer(final_answer_to_check)
print(result)