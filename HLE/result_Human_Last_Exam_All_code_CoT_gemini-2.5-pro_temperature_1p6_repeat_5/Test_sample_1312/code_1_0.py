import textwrap

def solve_genomics_question():
    """
    Analyzes a multiple-choice question about genome architecture and provides the correct answer with a detailed explanation.
    """
    question = "What intrinsic genome architectural feature is hypothesized to be a compensatory mechanism for preventing genetic deterioration in populations subjected to limited recombination?"
    
    options = {
        'A': "Tandem repeats",
        'B': "Chromosomal inversions",
        'C': "Transposable elements",
        'D': "Multigene families",
        'E': "Polyploidy"
    }
    
    correct_answer_key = 'B'
    
    explanation = textwrap.dedent(f"""
    The Problem:
    Populations with limited recombination (e.g., asexual or highly self-fertilizing species) face 'genetic deterioration'. This includes the accumulation of harmful mutations (Muller's Ratchet) and the negative effects of increasing homozygosity, which exposes recessive deleterious alleles (inbreeding depression).

    The Compensatory Mechanism:
    A chromosomal inversion is the reversal of a chromosome segment. In individuals heterozygous for an inversion, recombination is suppressed within the inverted region. This suppression has a critical consequence: it maintains heterozygosity for all the genes located within that large chromosomal block.

    Why This Helps:
    By maintaining heterozygosity, the inversion prevents harmful recessive alleles from being expressed. A single functional allele can mask its deleterious recessive counterpart on the other chromosome. This directly counteracts the genetic deterioration that would result from increased homozygosity, thus acting as a powerful compensatory mechanism.

    Why Other Options Are Less Likely:
    - Multigene families (D) and Polyploidy (E) provide compensation through redundancy, which is another valid mechanism, but the role of inversions in preserving heterozygosity is a distinct and major architectural adaptation.
    - Tandem repeats (A) and Transposable elements (C) are sources of mutation and are not considered primary compensatory mechanisms against genome-wide decay.
    """).strip()

    print("--- The Question ---")
    print(textwrap.fill(question, width=80))
    print("\n--- The Options ---")
    for key, value in options.items():
        print(f"{key}. {value}")
        
    print("\n--- The Correct Answer and Explanation ---")
    print(f"Correct Answer: {correct_answer_key}. {options[correct_answer_key]}")
    print("\nExplanation:")
    print(explanation)

# Execute the function to get the answer
solve_genomics_question()