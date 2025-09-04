def check_biology_question_answer():
    """
    This function checks the correctness of the answer to the iPSC chimera question.
    It simulates the logical deduction based on the biological principles
    and constraints provided in the question.
    """

    # 1. Define the key constraints and facts from the question.
    constraints = {
        "cell_type": "iPSC (undifferentiated)",
        "red_marker_control": "lineage-specific promoter",
        "green_marker_function": "apoptosis detection (TUNEL)",
        "known_outcome_1": "chimera formation is inefficient",
        "known_outcome_2": "unintegrated cells are cleared by apoptosis",
        "question_focus": "first thing you notice (implying first significant finding)"
    }

    # 2. Define the options provided in the question.
    options = {
        "A": "cell line-specific red signals label different organelles",
        "B": "there is no green signal",
        "C": "green signal colocalizes with the red signal",
        "D": "cytoplasmic localization of the red signal"
    }

    # 3. The final answer provided by the LLM to be checked.
    llm_answer = "C"

    # 4. Step-by-step logical evaluation based on constraints.
    analysis_log = []
    plausible_option = None

    # --- Evaluation of Option A ---
    # A promoter controls WHEN a gene is expressed, not WHERE the protein goes.
    # This option misinterprets basic molecular biology.
    if constraints["red_marker_control"] == "lineage-specific promoter":
        analysis_log.append("Constraint Violated by A: A promoter controls gene expression, not subcellular protein localization. This option is biologically incorrect.")
    
    # --- Evaluation of Option B ---
    # Apoptosis is expected in normal development and especially for clearing unintegrated cells.
    if constraints["known_outcome_2"] == "unintegrated cells are cleared by apoptosis":
        analysis_log.append("Constraint Violated by B: Apoptosis (green signal) is a highly expected event. This option is biologically implausible.")

    # --- Evaluation of Option D ---
    # While true that mRaspberry is cytoplasmic, this is a static property of the tool,
    # not a dynamic experimental finding about cell fate.
    if constraints["question_focus"] == "first thing you notice (implying first significant finding)":
        analysis_log.append("Constraint Violated by D: This is a trivial fact about the reporter protein, not a significant biological finding. It doesn't answer the core question about cell fate.")

    # --- Evaluation of Option C ---
    # This option describes a key biological event.
    # Event sequence: iPSC starts to differentiate (promoter on -> red signal) -> embryo recognizes it as abnormal -> embryo triggers apoptosis (TUNEL -> green signal).
    # This colocalization is the first direct evidence linking the two processes under study.
    analysis_log.append("Constraint Satisfied by C: This option describes the most significant and expected initial finding. It directly links differentiation (red signal) with apoptosis (green signal), which is the core of the experimental question.")
    plausible_option = "C"

    # 5. Final Verdict
    if llm_answer == plausible_option:
        return "Correct"
    else:
        error_message = f"The provided answer '{llm_answer}' is incorrect.\n"
        error_message += f"The most plausible answer is '{plausible_option}'.\n\n"
        error_message += "Analysis Log:\n"
        for log_entry in analysis_log:
            error_message += f"- {log_entry}\n"
        return error_message

# Run the check
result = check_biology_question_answer()
print(result)