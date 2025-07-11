def solve_genomics_question():
    """
    Analyzes the options for a genomics question and prints the correct answer and rationale.
    """
    question = "What intrinsic genome architectural feature is hypothesized to be a compensatory mechanism for preventing genetic deterioration in populations subjected to limited recombination?"

    options = {
        'A': 'Tandem repeats',
        'B': 'Chromosomal inversions',
        'C': 'Transposable elements',
        'D': 'Multigene families',
        'E': 'Polyploidy'
    }

    rationales = {
        'A': "Incorrect. Tandem repeats are sources of high mutation rates but do not provide a general compensatory mechanism against genetic deterioration.",
        'B': "Incorrect. Chromosomal inversions suppress recombination, exacerbating the problem rather than compensating for it.",
        'C': "Incorrect. Transposable elements can be mutagenic and their accumulation is often a sign of ineffective selection, not a compensatory mechanism.",
        'D': "Correct. Multigene families provide redundancy, so a mutation in one copy is not catastrophic. More importantly, they allow for gene conversion, a process that can repair mutated gene copies using other copies as templates, directly counteracting the accumulation of deleterious mutations (Muller's Ratchet).",
        'E': "Incorrect. While polyploidy provides redundancy, multigene families and the associated mechanism of gene conversion are considered a more direct and active architectural feature for compensating for the lack of recombination."
    }

    correct_answer_key = 'D'

    print(f"Question: {question}\n")
    print("Analyzing the options...")
    print("-------------------------")
    print(f"Chosen Answer: [{correct_answer_key}] {options[correct_answer_key]}")
    print(f"Rationale: {rationales[correct_answer_key]}")
    print("-------------------------")

solve_genomics_question()