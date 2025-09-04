import re

def check_answer_correctness(llm_response: str) -> str:
    """
    Checks the correctness of the LLM's answer based on biological principles.

    The function codifies the logic of dominant-negative mutations in dimeric proteins
    and evaluates the chosen option against these rules.
    """
    
    # Define the options as presented in the problem.
    # The final answer block has a consistent mapping for A, B, C, D.
    options = {
        'A': "protein aggregation and loss-of-function phenotype",
        'B': "protein degradation and loss-of-function of the wild-type allele",
        'C': "change of protein conformation and gain-of-function phenotype",
        'D': "loss of protein dimerization and wild-type phenotype"
    }

    # Extract the final answer choice (e.g., 'B') from the response
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Error: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."
    
    chosen_answer_key = match.group(1)
    chosen_answer_text = options.get(chosen_answer_key)

    # --- Rule-based Evaluation ---

    # Rule 1: Must be a Loss-of-Function (LoF) phenotype.
    if "gain-of-function" in chosen_answer_text:
        return (f"Incorrect. The answer '{chosen_answer_key}' is wrong because a dominant-negative mutation "
                f"is a type of loss-of-function, not a gain-of-function.")

    # Rule 2: Must be a mutant phenotype, not wild-type.
    if "wild-type phenotype" in chosen_answer_text:
        return (f"Incorrect. The answer '{chosen_answer_key}' is wrong because a dominant mutation "
                f"results in a mutant phenotype in a heterozygote, not a wild-type phenotype.")

    # Rule 3: The mechanism must allow for interference via dimerization.
    if "loss of protein dimerization" in chosen_answer_text:
        return (f"Incorrect. The answer '{chosen_answer_key}' is wrong because if the mutant protein completely "
                f"lost its ability to dimerize, it could not interfere with the wild-type protein. "
                f"This would describe a recessive mutation, not a dominant-negative one.")

    # Rule 4 & 5: Check for precision and find the best option.
    # The best option uses the most precise definition of a dominant-negative effect.
    best_option_key = None
    for key, text in options.items():
        # The most precise definition is "loss-of-function of the wild-type allele".
        if "loss-of-function of the wild-type allele" in text:
            best_option_key = key
            break
    
    if chosen_answer_key == best_option_key:
        return "Correct"
    else:
        # The chosen answer passed the basic checks but wasn't the most precise one.
        return (f"Incorrect. The chosen answer '{chosen_answer_key}' is suboptimal. "
                f"While '{options[chosen_answer_key]}' describes a plausible scenario, "
                f"the best answer is '{best_option_key}' ('{options[best_option_key]}'). "
                f"This is because it contains the phrase 'loss-of-function of the wild-type allele', "
                f"which is the precise, textbook definition of a dominant-negative effect.")

# The final answer provided by the LLM to be checked
llm_final_answer = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

### Step 1: Deconstruct the Core Biological Information

*   **Protein Function:** The protein is a transcription factor subunit that must form a dimer (a pair) to be active. Dimerization is an essential step for its function.
*   **Mutation X (Context):** This is a recessive loss-of-function mutation. The key implication is that in a heterozygote (one normal allele, one mutant X allele), the cell functions normally. This tells us that having only 50% of the normal amount of functional protein is sufficient for a wild-type phenotype. This concept is known as **haplosufficiency**.
*   **Mutation Y (The Focus):** This is a **dominant-negative** mutation located in the **dimerization domain**.
    *   "Dominant" means it causes a mutant phenotype even when a normal, wild-type (WT) allele is present.
    *   "Negative" means the mutant protein product not only loses its own function but also actively interferes with or "poisons" the function of the normal protein produced by the WT allele.

### Step 2: Analyze the Molecular Mechanism of the Dominant-Negative Effect

1.  Since the cell is haplosufficient (50% function is enough), the mutant Y protein must do more than just be non-functional. It must actively reduce the function of the WT protein to below this 50% threshold.
2.  The mutation is in the dimerization domain, which is the strongest clue about the mechanism of interference.
3.  For the mutant Y protein to interfere with the WT protein, it must be able to physically interact with it. If the mutation completely prevented dimerization, the Y protein would be inert, and the WT proteins would be free to form functional WT-WT dimers. Due to haplosufficiency, this would result in a normal phenotype, making the mutation recessive, which contradicts the problem statement.
4.  Therefore, the core mechanism is that the mutant Y subunit **must still be able to bind to a WT subunit**, forming a non-functional **heterodimer (WT-Y)**.
5.  This formation of non-functional "poisoned" complexes sequesters the functional WT subunits, preventing them from forming functional WT-WT dimers and causing a loss-of-function phenotype.

### Step 3: Evaluate the Answer Choices Based on the Mechanism

*   **A) protein aggregation and loss-of-function phenotype:** This is a plausible mechanism. A mutation in an interaction domain can lead to misfolding and cause proteins to clump together, sequestering the WT protein. However, the term "aggregation" can be general, and we should check if a more precise option exists.

*   **B) protein degradation and loss-of-function of the wild-type allele:** This option is highly accurate and specific. It has two components:
    *   **"loss-of-function of the wild-type allele"**: This is the precise, textbook definition of a dominant-negative effect. The mutant protein's action effectively nullifies the function of the protein produced by the normal, wild-type allele.
    *   **"protein degradation"**: This is a very common and specific molecular fate for the aberrant WT-Y heterodimers. The cell's quality control machinery (e.g., the ubiquitin-proteasome system) recognizes malformed or non-functional protein complexes and targets them for destruction. This degradation actively removes the WT protein from the cell (along with the mutant it's bound to), providing a direct mechanism for the "loss-of-function of the wild-type allele."

*   **C) change of protein conformation and gain-of-function phenotype:** This is incorrect. A dominant-negative mutation is a specific type of **loss-of-function** mutation, not a gain-of-function.

*   **D) loss of protein dimerization and wild-type phenotype:** This is incorrect. As explained in Step 2, the mutant protein must be able to dimerize to interfere. Furthermore, a dominant mutation results in a mutant phenotype, not a wild-type one.

### Step 4: Final Conclusion

Comparing the plausible options (A and B), option B is superior. It provides the most precise description of the genetic effect ("loss-of-function of the wild-type allele") and a very specific and common molecular mechanism ("protein degradation") that explains this effect. It is a more complete and accurate description of the molecular phenotype than the more general term "protein aggregation". The vast majority of the provided candidate answers also converge on this conclusion.

<<<B>>>
"""

# Run the check
result = check_answer_correctness(llm_final_answer)
print(result)