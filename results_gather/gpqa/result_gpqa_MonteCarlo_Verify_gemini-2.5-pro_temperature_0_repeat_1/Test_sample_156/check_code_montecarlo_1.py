def check_retrovirus_kit_design_answer(answer: str):
    """
    Checks the correctness of the answer for designing a retrovirus diagnostic kit.

    Args:
        answer: The selected option ('A', 'B', 'C', or 'D').

    Returns:
        A string indicating "Correct" or the reason for being incorrect.
    """
    # Define the core scientific principles and constraints from the question.
    constraints = {
        "virus_type": "retrovirus",
        "genome_type": "RNA",
        "required_kit_type": "molecular",
        "goal": "quick and accurate diagnosis"
    }

    # Analyze each option based on the principles.
    analysis = {
        'A': {
            "correct": False,
            "reason": "A retrovirus has an RNA genome, not a DNA genome. Therefore, the first step cannot be 'DNA sequencing' of the virus. The crucial step of converting RNA to cDNA (reverse transcription) is missing."
        },
        'B': {
            "correct": False,
            "reason": "Identifying a virus based on symptoms is not specific enough to design a molecular diagnostic kit. Designing PCR primers requires precise genetic sequence information, which cannot be obtained from symptoms alone."
        },
        'C': {
            "correct": True,
            "reason": "This is the correct approach. For a retrovirus (RNA genome), 'cDNA sequencing' correctly implies reverse transcription of RNA to cDNA followed by sequencing. A 'real time PCR' (RT-qPCR) kit is a standard, rapid, and accurate molecular diagnostic tool."
        },
        'D': {
            "correct": False,
            "reason": "The question asks for a 'molecular diagnostic kit'. An ELISA kit detecting IgG antibodies is a serological/immunological test, not a molecular one, as it detects the host's immune response rather than the virus's genetic material."
        }
    }

    if answer not in analysis:
        return f"Invalid option '{answer}'. Please choose from A, B, C, or D."

    selected_analysis = analysis[answer]

    if selected_analysis["correct"]:
        return "Correct"
    else:
        return f"Incorrect. The answer '{answer}' is wrong because: {selected_analysis['reason']}"

# The given answer from the other LLM
llm_answer = "C"

# Run the checking code
result = check_retrovirus_kit_design_answer(llm_answer)
print(result)