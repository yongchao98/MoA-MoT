def check_answer(llm_answer):
    """
    Checks the correctness of the answer based on biological principles
    of yeast cell cycle and gene expression during shmoo formation.
    """

    # --- Model the Biological State ---
    # Shmoo formation implies G1 arrest and high transcriptional activity for mating genes.
    # The assay targets 'active chromatin', i.e., regions of active transcription.
    state = {
        "cell_cycle_phase": "G1_arrest",
        "dominant_process": "transcription",
        "inhibited_process": "replication"
    }

    # --- Model the Protein Complexes ---
    # We assign a relative abundance score based on the biological state.
    # A higher score means more abundant in the ChIP-MS of active chromatin.
    complexes = {
        "A": {"name": "pre-replication complex", "function": "replication"},
        "B": {"name": "pre-initiation complex", "function": "transcription"},
        "C": {"name": "nucleosome histone complex", "function": "chromatin_structure"},
        "D": {"name": "enhancer protein complex", "function": "transcription"}
    }

    abundance_scores = {}

    for key, props in complexes.items():
        # C: Nucleosomes are the most fundamental and ubiquitous part of any chromatin prep.
        # They will be the most abundant by far.
        if props["function"] == "chromatin_structure":
            abundance_scores[key] = 100
        # B & D: These are essential for the dominant process (transcription).
        # They will be highly abundant on active chromatin.
        elif props["function"] == state["dominant_process"]:
            abundance_scores[key] = 80
        # A: This complex is for the inhibited process (replication).
        # While present in G1, it's not part of the 'active' machinery driving the phenotype.
        # Thus, it will be the least abundant in an assay for active chromatin proteome.
        elif props["function"] == state["inhibited_process"]:
            abundance_scores[key] = 10
        else:
            abundance_scores[key] = 0 # Should not happen

    # --- Determine the Least Abundant Complex ---
    # Find the key (complex letter) with the minimum score.
    least_abundant_complex = min(abundance_scores, key=abundance_scores.get)

    # --- Validate the LLM's Answer ---
    if llm_answer == least_abundant_complex:
        return "Correct"
    else:
        expected = complexes[least_abundant_complex]['name']
        provided = complexes[llm_answer]['name']
        reason = (f"The answer '{llm_answer}' ({provided}) is incorrect. "
                  f"The question asks for the LEAST observed complex in active chromatin during G1 arrest. "
                  f"In this state, transcription is highly active, while DNA replication is inhibited. "
                  f"Therefore, the '{expected}' (A), which is involved in the inhibited process of replication, "
                  f"should be the least observed. The '{provided}' is directly involved in the highly active process of transcription "
                  f"(or is a core structural component like histones) and would be abundant.")
        return reason

# The answer provided by the other LLM
llm_answer_choice = "A"

# Run the check
result = check_answer(llm_answer_choice)
print(result)