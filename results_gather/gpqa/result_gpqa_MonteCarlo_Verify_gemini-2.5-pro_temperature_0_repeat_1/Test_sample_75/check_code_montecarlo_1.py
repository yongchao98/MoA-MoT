def check_biology_answer_correctness():
    """
    Checks the correctness of the LLM's answer by verifying it against the
    explicit constraints given in the biological question.
    """
    # --- Define Constraints from the Question ---
    constraints = {
        "cell_type": "iPSC",  # Pluripotent, not yet differentiated.
        "reporter_promoter": "lineage-specific",
        "observation_timing": "first thing you notice",  # Implies the earliest timepoint, 12h.
        "first_timepoint_hr": 12,
        "apoptosis_stain": "TUNEL-FITC",  # Green signal for apoptosis.
        "cell_label": "mRaspberry",  # Red signal for iPSCs.
    }

    llm_answer = "C" # The answer to be checked.
    # Option C: "green signal colocalizes with the red signal"

    # --- Logical Verification Step 1: Check for Red Signal Presence ---
    # The core of the problem lies in whether the red signal (mRaspberry) would be expressed
    # at the first time point (12 hours).
    # Principle: A lineage-specific promoter only becomes active after a cell
    # commits to and differentiates into that specific lineage.
    # Pluripotent iPSCs need time to receive signals, differentiate, and then transcribe
    # and translate a gene into a detectable amount of protein.
    
    is_red_signal_expected_at_12h = False
    
    reasoning = (
        f"The promoter for the red fluorescent protein is '{constraints['reporter_promoter']}'. "
        f"This type of promoter is inactive in undifferentiated, pluripotent '{constraints['cell_type']}'s. "
        f"The '{constraints['observation_timing']}' implies the earliest timepoint of {constraints['first_timepoint_hr']} hours. "
        "This is an insufficient amount of time for the injected iPSCs to integrate, differentiate into a specific lineage, "
        "and express a detectable level of the mRaspberry protein."
    )

    # --- Logical Verification Step 2: Evaluate Answer 'C' ---
    # Answer 'C' requires the presence of BOTH red and green signals, and for them to be in the same location.
    if not is_red_signal_expected_at_12h:
        return (
            "Incorrect. The provided answer 'C' is not the 'first thing you notice'.\n"
            f"Reason: {reasoning}\n"
            "Because no red signal is expected at the 12-hour mark, it is impossible to observe its colocalization with the green signal (from apoptosis). "
            "Therefore, answer 'C' violates the temporal constraint of the question. While this colocalization might be a key finding later in the 48-hour experiment, it cannot be the first observation."
        )

    # --- Logical Verification Step 3: Evaluate other options (for completeness) ---
    # This part of the logic is technically unreachable if the above check is correct,
    # but it demonstrates a full analysis.
    # If we incorrectly assume the red signal IS present, is 'C' the best answer?
    # - A) "label different organelles": Unlikely, mRaspberry is cytoplasmic unless tagged. Not the main finding.
    # - B) "no green signal": Incorrect. Apoptosis is expected in normal development and for non-integrated cells.
    # - D) "cytoplasmic localization": True, but a trivial observation about the tool, not a key experimental result.
    # - C) "green signal colocalizes with the red signal": This describes the fate of the injected cells, which is the central goal of the experiment.
    # So, if the timing constraint is ignored, 'C' would be the most significant finding. However, the constraint cannot be ignored.

    return "Correct"

# Execute the check
result = check_biology_answer_correctness()
print(result)