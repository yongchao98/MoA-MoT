import re

def check_correctness_of_sars_cov_2_answer(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer to the SARS-CoV-2 molecular biology question.

    The function verifies the factual accuracy of each statement and then checks if the
    LLM correctly identified an incorrect statement. The question is flawed as both
    statements A and C contain factual errors. A valid answer is one that identifies
    either A or C as incorrect.

    Args:
        question_text (str): The full text of the question including all statements.
        llm_answer_text (str): The full text of the LLM's answer, including the final <<<X>>> tag.

    Returns:
        str: "Correct" if the answer is valid, or a string explaining the error.
    """

    # Extract the final answer choice (e.g., 'A') from the llm_answer_text string
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure to parse the answer. The answer format <<<X>>> was not found."
    
    final_answer_choice = match.group(1)

    # --- Scientific Facts Database ---
    # This section encodes established scientific knowledge to check against the statements.

    facts = {
        'A': {
            'is_correct': False,
            'reason': "The nsp10/nsp14-ExoN complex is an exoribonuclease. Its function is to CLEAVE or DEGRADE RNA for proofreading. The statement claims it 'prevents the breakdown', which is the exact opposite of its known biochemical function."
        },
        'B': {
            'is_correct': True,
            'reason': "This is a standard and accurate description of the -1 programmed ribosomal frameshifting (-1 PRF) mechanism in coronaviruses, including the conservation of the signal between SARS-CoV and SARS-CoV-2."
        },
        'C': {
            'is_correct': False,
            'reason': "This statement contains two factual errors. First, the relationship between frameshifting rate and pseudoknot dynamics is not a simple 'linear correlation'. Second, single-molecule studies show the SARS-CoV-2 pseudoknot has a three-state unfolding pathway, not a two-state one as claimed for both viruses."
        },
        'D': {
            'is_correct': True,
            'reason': "This statement accurately describes the initial trigger of apoptosis by ORF3a via the extrinsic pathway (caspase-8 activation), as supported by scientific literature. While an oversimplification of the entire process, it is not factually incorrect."
        }
    }

    # --- Evaluation Logic ---

    # Identify all statements that are factually incorrect based on our database.
    incorrect_statements = {key for key, value in facts.items() if not value['is_correct']}

    # Check if the LLM's chosen answer is one of the incorrect statements.
    if final_answer_choice in incorrect_statements:
        # The LLM correctly identified a factually incorrect statement.
        # Since the question is ambiguous (having multiple incorrect options), any correct identification is a valid answer.
        return "Correct"
    else:
        # The LLM chose a statement that is actually correct.
        error_reasons = [f"Statement {key}: {value['reason']}" for key, value in facts.items() if not value['is_correct']]
        return (f"Incorrect. The provided answer is '{final_answer_choice}', but this statement is factually correct. "
                f"The question asks for the incorrect statement. The actual incorrect statements are {sorted(list(incorrect_statements))} "
                f"for the following reasons: {' '.join(error_reasons)}")


# The user's provided question and answer
question = """
All the following statements about the molecular biology of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS‑CoV‑2) are correct except

A) SARS-CoV-2 nsp10/nsp14-ExoN operates as heterodimers in a mismatch repair mechanism. The N-terminal ExoN domain of nsp14 could bind to nsp10 making an active exonuclease complex that prevents the breakdown of dsRNA.
B) Programmed ribosomal frameshifting creates two polyproteins near to 5` end of the genome by moving back by 1 nucleotide with the help of slippery nucleotides, and pseudoknot. The SARS-CoV-2 programmed ribosomal frameshifting mostly has the same conformation as the SARS-CoV programmed ribosomal frameshifting.
C) The rate of frameshifting in vitro is linearly correlated with the number of conformations that a pseudoknot can adopt. Both SARS-CoV and SARS-CoV-2 Programmed -1 Frameshift Signals show two conformations when under tension, similar to other pseudoknots that induce comparable frameshifting rates.
D) SARS-CoV-2 ORF3a has the ability to trigger caspase-8 activation/cleavage, without affecting the expression levels of Bcl-2. Caspase-8 activation is recognized as a characteristic feature of the extrinsic apoptotic pathway via death receptors, while Bcl-2 plays a crucial role in initiating the mitochondrial pathway. This suggests that the mechanism through which SARS-CoV-2 ORF3a induces apoptosis is via the extrinsic apoptotic pathway.
"""

llm_answer = """
Here is a step-by-step analysis to decide the final answer.

**1. Deconstruct and Analyze Each Statement**

The goal is to find the single statement that is factually incorrect.

*   **Statement A: SARS-CoV-2 nsp10/nsp14-ExoN operates as heterodimers in a mismatch repair mechanism. The N-terminal ExoN domain of nsp14 could bind to nsp10 making an active exonuclease complex that prevents the breakdown of dsRNA.**
    *   **Analysis:** The first part of the statement is correct. The nsp10/nsp14 complex is a heterodimer that acts as the virus's proofreading machine (mismatch repair) via the exoribonuclease (ExoN) activity of nsp14. The critical error is in the final clause: "...that **prevents** the breakdown of dsRNA."
    *   **Factual Check:** An exonuclease, by its very definition, is an enzyme that *cleaves* or *breaks down* nucleic acids. The function of this complex is to remove misincorporated nucleotides by actively degrading the RNA strand. Stating that it *prevents* breakdown is a direct contradiction of its fundamental biochemical function.
    *   **Conclusion:** This statement contains a clear and fundamental factual error.

*   **Statement B: Programmed ribosomal frameshifting creates two polyproteins near to 5` end of the genome by moving back by 1 nucleotide with the help of slippery nucleotides, and pseudoknot. The SARS-CoV-2 programmed ribosomal frameshifting mostly has the same conformation as the SARS-CoV programmed ribosomal frameshifting.**
    *   **Analysis:** This statement provides a standard and accurate description of the -1 programmed ribosomal frameshifting (-1 PRF) mechanism in coronaviruses. It correctly identifies the key elements (slippery sequence, pseudoknot) and the outcome (production of pp1a and pp1ab polyproteins). It also correctly notes the high degree of structural conservation of the frameshifting signal between SARS-CoV and SARS-CoV-2.
    *   **Conclusion:** This statement is factually correct.

*   **Statement C: The rate of frameshifting in vitro is linearly correlated with the number of conformations that a pseudoknot can adopt. Both SARS-CoV and SARS-CoV-2 Programmed -1 Frameshift Signals show two conformations when under tension, similar to other pseudoknots that induce comparable frameshifting rates.**
    *   **Analysis:** This statement makes two specific biophysical claims that are incorrect.
        1.  The relationship between frameshifting efficiency and pseudoknot dynamics is complex and not a simple "linear correlation."
        2.  More importantly, single-molecule force spectroscopy studies have revealed a key difference between the two viruses: the SARS-CoV pseudoknot unfolds in a two-state manner, but the SARS-CoV-2 pseudoknot unfolds via a **three-state pathway** that includes a stable intermediate conformation. The claim that *both* show two conformations is factually wrong.
    *   **Conclusion:** This statement contains at least two factual errors.

*   **Statement D: SARS-CoV-2 ORF3a has the ability to trigger caspase-8 activation/cleavage, without affecting the expression levels of Bcl-2... This suggests that the mechanism through which SARS-CoV-2 ORF3a induces apoptosis is via the extrinsic apoptotic pathway.**
    *   **Analysis:** The premises of this statement are largely supported by research. ORF3a does activate caspase-8 (a hallmark of the extrinsic pathway), and some studies show this occurs without altering Bcl-2 levels. The conclusion that this "suggests" an extrinsic pathway is a reasonable interpretation of the initial trigger. While the full process is more complex and involves crosstalk with the intrinsic pathway, the statement accurately describes the initiation events.
    *   **Conclusion:** This statement is considered factually correct, albeit an oversimplification of the entire process.

**2. Final Decision: Comparing the Incorrect Statements (A vs. C)**

Both statements A and C are factually incorrect. We must choose the one that is the intended "exception."

*   **Error in Statement A:** This is a fundamental contradiction of a core biochemical concept. It describes an enzyme (an exonuclease) as performing the exact opposite of its defined function (preventing breakdown vs. causing breakdown). This is an unambiguous and significant error.
*   **Error in Statement C:** This involves misstating specific, nuanced results from biophysical experiments (a three-state pathway is described as two-state) and oversimplifying a complex relationship (describing it as "linear"). While factually incorrect, it is an error in the details of a model rather than a reversal of a fundamental definition.

The error in Statement A is more foundational and represents a greater misunderstanding of molecular biology principles. An "exonuclease complex that prevents breakdown" is an oxymoron. Therefore, it is the most definitively incorrect statement.

<<<A>>>
"""

# Run the check
result = check_correctness_of_sars_cov_2_answer(question, llm_answer)
print(result)