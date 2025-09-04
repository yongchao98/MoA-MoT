def check_correctness():
    """
    This function checks the correctness of the provided answer for designing a molecular diagnostic kit for a retrovirus.

    The function evaluates each option based on key scientific principles:
    1.  **Kit Type**: The question asks for a "molecular" kit, which detects nucleic acids (DNA/RNA).
    2.  **Pathogen Biology**: The target is a "retrovirus," which has an RNA genome. Any DNA-based analysis (like PCR or sequencing) requires a preliminary step of reverse transcription to create complementary DNA (cDNA).
    3.  **Diagnostic Goal**: The kit must be "quick and accurate," which is crucial for managing an outbreak.
    """
    # The final answer provided by the LLM that we need to check.
    llm_answer = "A"

    # A dictionary to hold the analysis of each option.
    analysis = {
        "A": {
            "is_correct": False,
            "reason": ""
        },
        "B": {
            "is_correct": False,
            "reason": ""
        },
        "C": {
            "is_correct": False,
            "reason": ""
        },
        "D": {
            "is_correct": False,
            "reason": ""
        }
    }

    # --- Evaluate Option A ---
    # "First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis."
    # - Correctly identifies the need for cDNA from the RNA genome (satisfies pathogen biology).
    # - Proposes a molecular method (real-time PCR) that is quick and accurate (satisfies kit type and goal).
    analysis["A"]["is_correct"] = True
    analysis["A"]["reason"] = "This option correctly follows the standard scientific workflow: convert the retroviral RNA to cDNA for sequencing, then use that information to design a rapid and accurate molecular test (real-time PCR)."

    # --- Evaluate Option B ---
    # "First you identify the virus using the information about the symptoms... You then design a nested PCR kit..."
    # - Flaw: Symptoms are non-specific and cannot be used to determine a virus's genetic sequence, which is required to design PCR primers.
    analysis["B"]["is_correct"] = False
    analysis["B"]["reason"] = "The initial step of identifying a virus based on symptoms is not scientifically valid for designing a specific molecular test."

    # --- Evaluate Option C ---
    # "First you identify the virus by performing DNA sequencing. You then develop a PCR kit..."
    # - Flaw: A retrovirus has an RNA genome. Direct DNA sequencing of the viral particle is biologically incorrect and would not work. It misses the essential reverse transcription step.
    analysis["C"]["is_correct"] = False
    analysis["C"]["reason"] = "This option is incorrect because a retrovirus has an RNA genome, not a DNA genome. The first step must involve converting RNA to cDNA."

    # --- Evaluate Option D ---
    # "First you identify the IgG antibodies for the infection. You then develop an ELISA kit..."
    # - Flaw 1: This describes a serological test (detecting antibodies), not a molecular test (detecting genetic material).
    # - Flaw 2: IgG antibodies appear late in an infection, making an IgG-based test unsuitable for "quick" diagnosis during an outbreak's early stages.
    analysis["D"]["is_correct"] = False
    analysis["D"]["reason"] = "This is not a molecular kit as requested. Furthermore, detecting IgG antibodies is not suitable for quick, early diagnosis of an active infection."

    # --- Final Verdict ---
    correct_option = None
    for option, details in analysis.items():
        if details["is_correct"]:
            correct_option = option
            break

    if llm_answer == correct_option:
        return "Correct"
    else:
        return f"Incorrect. The provided answer '{llm_answer}' is wrong. The correct answer is '{correct_option}' because {analysis[correct_option]['reason']}. The reason answer '{llm_answer}' is wrong is: {analysis[llm_answer]['reason']}"

# Execute the check and print the result
print(check_correctness())