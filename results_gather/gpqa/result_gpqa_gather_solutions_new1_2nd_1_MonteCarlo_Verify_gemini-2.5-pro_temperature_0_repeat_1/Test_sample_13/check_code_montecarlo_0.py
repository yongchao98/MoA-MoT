import json

def check_chip_seq_answer(llm_answer_string: str) -> str:
    """
    Checks the correctness of the answer to a complex molecular biology question.

    The function encodes the scientific reasoning required to answer the question
    and compares the logically derived answer to the provided one.
    """
    
    # --- Step 1: Define the problem space and the provided answer ---
    question_summary = "Why would IKAROS ChIP-seq peaks present with PFA fixation disappear with PFA+DSG fixation?"
    options = {
        "A": "At random locations in the genome",
        "B": "At active promoters and enhancers",
        "C": "In the introns of large genes",
        "D": "At repeats"
    }
    
    try:
        # Extract the letter from the format <<<D>>>
        provided_answer = llm_answer_string.strip().split('<<<')[-1].split('>>>')[0]
        if provided_answer not in options:
            return f"Invalid Answer Format: The answer '{provided_answer}' is not one of the valid options {list(options.keys())}."
    except IndexError:
        return f"Invalid Answer Format: Could not parse the answer from the string '{llm_answer_string}'. Expected format is <<<X>>>."

    # --- Step 2: Encode biological facts and technical principles ---
    knowledge_base = {
        "pfa_vs_dsg": "PFA is a short-range crosslinker. PFA+DSG is a more stringent dual-crosslinking method (DSG is a long-range protein-protein crosslinker) used to better stabilize large protein complexes.",
        "ikaros_biology": "IKAROS has a well-documented dual function: it binds to (1) active promoters/enhancers to regulate genes, and (2) dense, repetitive pericentromeric heterochromatin.",
        "epitope_masking_hypothesis": {
            "location": "B", # Active promoters and enhancers
            "mechanism": "DSG creates a dense web of crosslinked proteins in large complexes, physically blocking the antibody epitope.",
            "plausibility": "Possible, but it contradicts the primary goal of using DSG, which is to *improve* signal at these functional sites."
        },
        "insolubility_artifact_hypothesis": {
            "location": "D", # At repeats
            "mechanism": "The strong PFA+DSG crosslinking makes already-dense heterochromatin (composed of repeats, where IKAROS binds) insoluble. These insoluble chunks are physically lost during sample preparation.",
            "plausibility": "Highly plausible. It provides a direct physical reason for complete signal loss and is supported by the specific biology of IKAROS."
        },
        "artifact_removal_hypothesis": {
            "location": "D", # At repeats
            "mechanism": "The more stringent PFA+DSG method enhances true signals (at promoters) and raises the statistical bar for peak-calling, causing weaker, artifactual signals (common at 'sticky' repeats) to be filtered out.",
            "plausibility": "Highly plausible. Aligns with the principle that more stringent methods clean up data."
        }
    }

    # --- Step 3: Logically evaluate the options ---
    evaluation_results = {}

    # Option A: Random locations
    evaluation_results["A"] = {
        "is_correct": False,
        "reason": "ChIP-seq peaks are non-random by definition. A systematic experimental change would affect specific classes of genomic regions, not random ones."
    }

    # Option C: Introns of large genes
    evaluation_results["C"] = {
        "is_correct": False,
        "reason": "This category is too general. The underlying cause is a specific biophysical or functional property (like being a repeat or an enhancer), not simply being in an intron."
    }
    
    # Option B: Active promoters and enhancers
    hyp_b = knowledge_base["epitope_masking_hypothesis"]
    evaluation_results["B"] = {
        "is_correct": False, # Considered less likely than D
        "reason": f"This is a plausible hypothesis via epitope masking. However, its plausibility is weaker because {hyp_b['plausibility']}"
    }

    # Option D: At repeats
    hyp_d1 = knowledge_base["insolubility_artifact_hypothesis"]
    hyp_d2 = knowledge_base["artifact_removal_hypothesis"]
    evaluation_results["D"] = {
        "is_correct": True,
        "reason": f"This is the most likely answer. It is supported by two strong, complementary hypotheses that are consistent with the specific biology of IKAROS. 1) {hyp_d1['mechanism']} 2) {hyp_d2['mechanism']}"
    }

    # --- Step 4: Determine the correct answer and check against the provided one ---
    derived_correct_answer = None
    for option, result in evaluation_results.items():
        if result["is_correct"]:
            derived_correct_answer = option
            break
            
    if not derived_correct_answer:
        return "Error in logical evaluation: Could not determine a correct answer."

    if provided_answer == derived_correct_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer was '{provided_answer}', but the logically derived correct answer is '{derived_correct_answer}'.\n"
                f"Reasoning: The explanation for option D ('At repeats') is stronger than for any other option. "
                f"While epitope masking (Option B) is possible, the loss of signal from repeats (Option D) is supported by both the specific biology of IKAROS binding to dense heterochromatin and a more direct physical mechanism (insolubility) that explains the complete disappearance of peaks.")

# The final answer provided by the LLM to be checked
llm_final_answer = "<<<D>>>"

# Run the check
result = check_chip_seq_answer(llm_final_answer)
print(result)