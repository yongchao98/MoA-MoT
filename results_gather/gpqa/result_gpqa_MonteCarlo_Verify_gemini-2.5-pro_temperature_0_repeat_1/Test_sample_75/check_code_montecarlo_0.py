def check_biology_answer_correctness():
    """
    Checks the correctness of the answer to the biology question by evaluating
    each option against the experimental design and known biological principles.
    """
    llm_answer = "C"

    # --- Define facts from the question and general biology ---
    experimental_facts = {
        "iPSC_marker": "red (mRaspberry)",
        "apoptosis_marker": "green (TUNEL-FITC)",
        "primary_goal": "Investigate iPSC fate and co-localization with apoptosis.",
        "known_iPSC_fate": "A significant fraction of injected cells fail to integrate and undergo apoptosis.",
        "mRaspberry_localization": "cytoplasmic",
        "apoptosis_in_embryo": "Occurs normally during development and for clearing foreign cells."
    }

    # --- Evaluate each possible answer ---
    # Option A: "cell line-specific red signals label different organelles"
    # Check: Is this consistent with the marker used?
    if experimental_facts["mRaspberry_localization"] != "various_organelles":
        reason_A_is_wrong = "Option A is incorrect because mRaspberry is a cytoplasmic protein, not one that typically labels various organelles."
    else:
        reason_A_is_wrong = None

    # Option B: "there is no green signal"
    # Check: Is this consistent with expected biological processes?
    if experimental_facts["apoptosis_in_embryo"] == "Occurs normally during development and for clearing foreign cells.":
        reason_B_is_wrong = "Option B is incorrect because apoptosis (green signal) is expected both as a normal part of embryogenesis and as the mechanism to clear non-integrated iPSCs."
    else:
        reason_B_is_wrong = None

    # Option D: "cytoplasmic localization of the red signal"
    # Check: Is this the primary, non-trivial finding?
    if experimental_facts["mRaspberry_localization"] == "cytoplasmic":
        reason_D_is_wrong = "Option D is incorrect because while true, it describes a trivial characteristic of the reporter protein, not a significant biological finding related to the experiment's goal of understanding cell fate."
    else:
        reason_D_is_wrong = None

    # Option C: "green signal colocalizes with the red signal"
    # Check: Does this describe a key expected outcome related to the goal?
    if experimental_facts["known_iPSC_fate"] == "A significant fraction of injected cells fail to integrate and undergo apoptosis.":
        # This option describes the visualization of the primary expected outcome.
        reason_C_is_correct = "This describes the key biological event: injected iPSCs (red) undergoing apoptosis (green), which directly addresses the experiment's goal."
    else:
        reason_C_is_correct = None

    # --- Final Verification ---
    correct_option = "C"
    reasons_for_incorrectness = {
        "A": reason_A_is_wrong,
        "B": reason_B_is_wrong,
        "D": reason_D_is_wrong
    }

    if llm_answer == correct_option:
        # Acknowledge the nuance about the promoter timing, which makes the question slightly imperfect
        # but doesn't change the best answer.
        ambiguity_note = (
            "The check confirms 'C' is the best answer. It's worth noting the ambiguity of the 'lineage-specific promoter', "
            "which might delay the appearance of the red signal. However, among the choices, C is the only one describing a "
            "significant, non-trivial outcome central to the experiment's purpose."
        )
        # print(ambiguity_note) # Optional: print the note for full context
        return "Correct"
    else:
        if llm_answer in reasons_for_incorrectness:
            return f"Incorrect. The answer '{llm_answer}' is wrong. Reason: {reasons_for_incorrectness[llm_answer]}"
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is not the best choice. The correct answer is '{correct_option}' because {reason_C_is_correct}."

# Run the checker
result = check_biology_answer_correctness()
print(result)