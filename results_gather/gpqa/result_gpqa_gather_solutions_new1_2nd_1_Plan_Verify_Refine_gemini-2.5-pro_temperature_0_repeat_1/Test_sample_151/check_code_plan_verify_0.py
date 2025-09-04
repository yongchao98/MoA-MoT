import json

def check_answer():
    """
    Checks the correctness of the LLM's answer to the biology question.
    """
    question_options = {
        "A": "enhancer protein complex",
        "B": "pre-replication complex",
        "C": "pre-initiation complex",
        "D": "nucleosome histone complex"
    }

    # The provided answer from the LLM
    llm_answer_letter = "B"

    # Step 1: Define the biological state based on the question.
    # Shmoo formation in yeast means G1 arrest.
    cellular_state = {
        "cell_cycle_phase": "G1 arrest",
        "transcription": "highly active",
        "dna_replication": "inhibited"
    }

    # Step 2: Define the function and expected abundance of each complex in this state.
    # The experiment targets "active chromatin", so abundance is relative to this activity.
    complex_evaluation = {
        "enhancer protein complex": {
            "function": "Drives active transcription of specific genes.",
            "expected_abundance": "High",
            "reason": "The cell is undergoing massive transcriptional reprogramming for mating."
        },
        "pre-replication complex": {
            "function": "Licenses DNA for replication, a process that occurs in S phase.",
            "expected_abundance": "Low",
            "reason": "The cell is arrested in G1 and DNA replication is actively inhibited. This complex is associated with a stalled process."
        },
        "pre-initiation complex": {
            "function": "Includes RNA polymerase; essential for starting transcription.",
            "expected_abundance": "High",
            "reason": "The cell is highly transcriptionally active."
        },
        "nucleosome histone complex": {
            "function": "Fundamental structural unit of all chromatin (active and inactive).",
            "expected_abundance": "Very High",
            "reason": "It is the most ubiquitous protein complex on DNA."
        }
    }

    # Step 3: Determine the correct answer based on the evaluation.
    # The question asks for the LEAST observed complex.
    least_abundant_complex = None
    min_abundance_level = float('inf')
    abundance_scale = {"Low": 1, "High": 2, "Very High": 3}

    for complex_name, details in complex_evaluation.items():
        level = abundance_scale[details["expected_abundance"]]
        if level < min_abundance_level:
            min_abundance_level = level
            least_abundant_complex = complex_name

    # Find the letter corresponding to the correct complex
    correct_letter = None
    for letter, name in question_options.items():
        if name == least_abundant_complex:
            correct_letter = letter
            break

    # Step 4: Compare the derived correct answer with the LLM's answer.
    if llm_answer_letter == correct_letter:
        return "Correct"
    else:
        reasoning = {
            "llm_answer_letter": llm_answer_letter,
            "llm_answer_complex": question_options.get(llm_answer_letter, "Invalid Option"),
            "correct_letter": correct_letter,
            "correct_complex": least_abundant_complex,
            "explanation": f"The question asks for the least observed complex in a G1-arrested, transcriptionally active cell. The process of DNA replication is inhibited in this state. Therefore, the '{least_abundant_complex}', which is responsible for preparing for replication, would be the least observed. This corresponds to option {correct_letter}. The provided answer was {llm_answer_letter} ('{question_options.get(llm_answer_letter, 'Invalid Option')}'), which is incorrect."
        }
        return f"Incorrect. Reason: {json.dumps(reasoning, indent=2)}"

# Run the check
result = check_answer()
print(result)