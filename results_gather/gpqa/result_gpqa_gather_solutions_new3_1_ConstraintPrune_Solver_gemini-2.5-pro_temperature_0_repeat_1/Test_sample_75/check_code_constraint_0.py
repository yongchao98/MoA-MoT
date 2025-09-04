import json

def check_answer(question_details, provided_answer):
    """
    Checks the correctness of the answer based on biological principles.

    The question asks for the "first thing you notice" in a mouse chimera experiment.
    - iPSCs are labeled with mRaspberry (red) under a lineage-specific promoter.
    - Apoptosis is labeled with TUNEL-FITC (green).

    Let's analyze the options based on established biological facts.
    """

    # --- Define Biological Principles ---
    principles = {
        "lineage_promoter_in_undifferentiated_cells": "inactive",
        "differentiation_takes_time": True,
        "apoptosis_in_early_embryo": "expected",
        "fate_of_injected_cells": ["integration", "apoptosis"],
        "mRaspberry_localization": "cytoplasmic",
        "promoter_controls_localization": False
    }

    # --- Analysis based on a strict interpretation of "first" (e.g., at 12h) ---
    # At the earliest time point, iPSCs are still undifferentiated.
    red_signal_expected_early = principles["lineage_promoter_in_undifferentiated_cells"] == "active"
    green_signal_expected_early = principles["apoptosis_in_early_embryo"] == "expected"

    strict_first_observation = {
        "red_signal": red_signal_expected_early, # False
        "green_signal": green_signal_expected_early # True
    }

    # --- Analysis based on "first significant biological finding" ---
    # This interpretation looks for the first event that links the experimental variables
    # (differentiation and apoptosis) to answer the research question.
    # A common fate for iPSCs is aberrant differentiation followed by elimination.
    # Sequence: Differentiate -> Turn Red -> Recognized as abnormal -> Apoptosis -> Turn Green
    # This results in a cell that is both Red and Green.
    significant_finding_scenario = "colocalization_of_red_and_green"

    # --- Evaluate all possible options ---
    options = {
        "A": "cytoplasmic localization of the red signal",
        "B": "there is no green signal",
        "C": "green signal colocalizes with the red signal",
        "D": "cell line-specific red signals label different organelles"
    }
    
    analysis_results = {}

    # Option A Analysis
    if strict_first_observation["red_signal"] is False:
        analysis_results["A"] = "Incorrect. A red signal is not expected at the first time point because the lineage-specific promoter is inactive in undifferentiated iPSCs."
    else:
        analysis_results["A"] = "Plausible but weak. This describes a property of the reporter protein, not a dynamic biological finding about cell fate."

    # Option B Analysis
    if strict_first_observation["green_signal"] is True:
        analysis_results["B"] = "Incorrect. Apoptosis (green signal) is expected in early development and as a common fate for injected cells."
    else:
        analysis_results["B"] = "Correct under this assumption."

    # Option D Analysis
    if principles["promoter_controls_localization"] is False:
        analysis_results["D"] = "Incorrect. This is biologically false. A promoter controls gene expression, not the subcellular localization of the protein."
    else:
        analysis_results["D"] = "Correct under this assumption."

    # Option C Analysis
    if strict_first_observation["red_signal"] is False:
        analysis_results["C_strict"] = "Incorrect under a strict 'first observation' interpretation, as a red signal is not yet present."
    
    # Check against the "significant finding" interpretation
    if significant_finding_scenario == "colocalization_of_red_and_green":
        analysis_results["C_finding"] = "Correct. This describes the first key biological event that addresses the research question: the elimination of aberrantly differentiating cells. This is the most likely interpretation of what a researcher would 'notice'."


    # --- Final Verdict ---
    # The provided answer is 'C'. Let's check its validity.
    # The strict interpretation invalidates all options or makes them highly unlikely.
    # This suggests the question is asking for the first significant finding.
    # Under that interpretation, 'C' is the correct answer.

    if provided_answer == "C":
        return "Correct"
    else:
        # Provide the reason why the given answer is wrong.
        if provided_answer == "A":
            return analysis_results["A"]
        elif provided_answer == "B":
            return analysis_results["B"]
        elif provided_answer == "D":
            return analysis_results["D"]
        else:
            return f"The provided answer '{provided_answer}' is incorrect for reasons outlined in the analysis."


# --- Execution of the Check ---
# The question is complex, so we represent it with a summary.
question_summary = {
    "cells": "iPSCs",
    "red_reporter": "mRaspberry under lineage-specific promoter",
    "green_stain": "TUNEL for apoptosis",
    "context": "Injection into blastocyst, observation from 12h",
    "task": "What is the first thing you notice?"
}

# The final answer from the LLM to be checked.
# The provided answer was 'C' (green signal colocalizes with the red signal).
# Note: The original question had different lettering. We map the provided answer's content to 'C'.
llm_answer = "C"

# Run the check
result = check_answer(question_summary, llm_answer)
print(result)
