def check_answer_correctness():
    """
    This function checks the correctness of the provided answer for designing a
    diagnostic kit for a retrovirus.
    """
    question_constraints = {
        "organism_type": "retrovirus",
        "genome_type": "RNA",
        "kit_type": "molecular",
        "speed": "quick"
    }

    llm_answer = "B"

    # Analysis of each option based on scientific principles
    options_analysis = {
        "A": {
            "correct": False,
            "reason": "This describes a serological test (ELISA for IgG antibodies), not a molecular diagnostic kit. A molecular kit detects the pathogen's genetic material. Furthermore, IgG antibodies appear later in an infection, making this method unsuitable for 'quick detection' in a new outbreak."
        },
        "B": {
            "correct": True,
            "reason": "This is the correct approach. A retrovirus has an RNA genome, so it must be converted to cDNA (complementary DNA) via reverse transcription before sequencing. Real-time PCR (RT-qPCR) is the gold standard for rapid, sensitive, and quantitative detection of viral RNA, making it perfect for a 'quick and accurate' molecular diagnostic kit."
        },
        "C": {
            "correct": False,
            "reason": "Identifying a virus based on symptoms is unreliable for designing a specific molecular test, as different viruses can cause similar symptoms. A definitive genetic sequence is required. While nested PCR is sensitive, it is generally slower and more prone to contamination than real-time PCR."
        },
        "D": {
            "correct": False,
            "reason": "This is fundamentally incorrect because a retrovirus has an RNA genome, not a DNA genome. The first step cannot be 'DNA sequencing'. The RNA must first be converted to cDNA."
        }
    }

    # Check if the LLM's answer is the correct one
    if llm_answer in options_analysis:
        if options_analysis[llm_answer]["correct"]:
            return "Correct"
        else:
            # Return the reason why the selected answer is wrong
            return options_analysis[llm_answer]["reason"]
    else:
        return f"Invalid answer option provided: {llm_answer}"

# Execute the check and print the result
result = check_answer_correctness()
print(result)