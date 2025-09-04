import collections

def check_sars_cov2_answer_correctness():
    """
    Verifies the correctness of the given answer by checking each statement
    against established scientific facts about SARS-CoV-2 molecular biology.
    """
    
    # The provided answer to the question.
    llm_answer = "D"

    # A database of verified facts for each statement.
    # Each entry contains a boolean indicating if the statement is correct
    # and a reason explaining why.
    FactCheckResult = collections.namedtuple('FactCheckResult', ['is_correct', 'reason'])
    facts = {
        "A": FactCheckResult(
            is_correct=True,
            reason="Statement A is correct. Single-molecule studies (e.g., Neupane et al., 2021, Nature Communications) have shown that the SARS-CoV-2 frameshifting pseudoknot exhibits dynamic switching between at least two conformations under tension, and this behavior is directly related to frameshifting efficiency."
        ),
        "B": FactCheckResult(
            is_correct=True,
            reason="Statement B is correct. Research (e.g., Ren et al., 2020, Cell Death & Differentiation) has demonstrated that the SARS-CoV-2 ORF3a protein induces apoptosis by activating the extrinsic pathway, confirmed by caspase-8 cleavage, without altering the levels of the intrinsic pathway protein Bcl-2."
        ),
        "C": FactCheckResult(
            is_correct=True,
            reason="Statement C is correct. This is an accurate description of the canonical -1 programmed ribosomal frameshifting (-1 PRF) mechanism in coronaviruses, which occurs at the ORF1a/1b junction and requires both a slippery sequence and a downstream pseudoknot. The high structural conservation between SARS-CoV and SARS-CoV-2 frameshifting elements is also well-documented."
        ),
        "D": FactCheckResult(
            is_correct=False,
            reason="Statement D is incorrect because it misrepresents the function of the nsp10/nsp14 complex. This complex is a 3'-to-5' exoribonuclease (ExoN) that acts as a proofreader. Its function is to *degrade* and *remove* mismatched nucleotides from the newly synthesized RNA strand to correct errors. It *causes* RNA breakdown for fidelity, it does not *prevent* the breakdown of dsRNA."
        )
    }

    # The question asks to identify the INCORRECT statement.
    # We need to find which statement our fact-check identifies as incorrect.
    truly_incorrect_statement = None
    for statement_key, result in facts.items():
        if not result.is_correct:
            truly_incorrect_statement = statement_key
            break
            
    # Compare the LLM's answer with our verified incorrect statement.
    if llm_answer == truly_incorrect_statement:
        return "Correct"
    else:
        if truly_incorrect_statement is None:
            # This case should not happen given the facts.
            return f"The provided answer '{llm_answer}' is incorrect because all statements were found to be correct."
        else:
            # The LLM identified the wrong statement.
            reason_for_llm_error = facts[llm_answer].reason
            reason_for_correct_answer = facts[truly_incorrect_statement].reason
            return (f"The provided answer '{llm_answer}' is incorrect. "
                    f"Statement '{llm_answer}' is actually correct. {reason_for_llm_error}. "
                    f"The truly incorrect statement is '{truly_incorrect_statement}'. {reason_for_correct_answer}")

# Execute the check and print the result.
result = check_sars_cov2_answer_correctness()
print(result)