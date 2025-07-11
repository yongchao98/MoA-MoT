import textwrap

def analyze_gene_duplication_fates():
    """
    Analyzes the mechanisms responsible for the retention and divergence of
    duplicate genes and prints the most likely answer.
    """
    question = "Which of the following mechanisms is most likely responsible for the retention and divergence of duplicate genes in eukaryotic genomes?"

    # A dictionary to hold the answer choices and their detailed explanations.
    choices = {
        'A': {
            'term': "Gene conversion",
            'explanation': "A process that makes two duplicate genes more similar by copying a sequence from one to the other. It actively works AGAINST divergence."
        },
        'B': {
            'term': "Pseudogenization",
            'explanation': "The most common fate, where one gene copy accumulates mutations and becomes non-functional. While this involves divergence, it doesn't lead to the retention of two functional genes, but rather one functional gene and one pseudogene."
        },
        'C': {
            'term': "Neofunctionalization",
            'explanation': "One gene copy retains the original function, while the other copy accumulates new mutations and acquires a novel function. This provides a strong evolutionary reason to retain both copies and requires them to diverge."
        },
        'D': {
            'term': "Subfunctionalization",
            'explanation': "The ancestral gene had multiple functions, and after duplication, each gene copy specializes in performing a subset of these original functions. This also explains both retention (both copies are now required) and divergence (as they specialize)."
        },
        'E': {
            'term': "Adaptive radiation",
            'explanation': "A macroevolutionary pattern of rapid diversification of species, not a molecular mechanism for retaining individual genes within a genome."
        }
    }

    print("Analyzing the following question:\n")
    print(textwrap.fill(question, width=80))
    print("\n----------------------------------\n")
    print("Evaluating the potential mechanisms:\n")

    best_choices = []

    for key, value in choices.items():
        # Evaluate based on the core requirements: retention and divergence.
        promotes_retention = "Neofunctionalization" in value['term'] or "Subfunctionalization" in value['term']
        promotes_divergence = "Neofunctionalization" in value['term'] or "Subfunctionalization" in value['term'] or "Pseudogenization" in value['term']
        
        print(f"Choice {key}: {value['term']}")
        print(textwrap.fill(f"  - Explanation: {value['explanation']}", width=80, initial_indent="  ", subsequent_indent="  "))
        
        if promotes_retention and promotes_divergence:
            print("  - Verdict: This is a very strong candidate as it explains both retention and divergence of functional genes.")
            best_choices.append(key)
        elif "Pseudogenization" in value['term']:
            print("  - Verdict: A very common fate, but it results in the loss of one functional copy, which doesn't fully fit the 'retention of duplicate genes' framing.")
        else:
            print("  - Verdict: This is not a primary mechanism for the retention AND divergence of duplicate genes.")
        print("-" * 30)

    print("\nConclusion:\n")
    conclusion_text = (
        "Both Neofunctionalization (C) and Subfunctionalization (D) are the primary accepted models explaining the "
        "retention and divergence of duplicate functional genes. Pseudogenization (B) is the most common fate but involves "
        "the loss of function. Gene conversion (A) opposes divergence, and Adaptive Radiation (E) is a process at the species level, not the gene level.\n\n"
        "Between Neofunctionalization and Subfunctionalization, both are highly plausible. However, they are often presented as the two key mechanisms responsible for this phenomenon. "
        "Both are considered 'most likely' depending on the specific context. In the context of a multiple-choice question, "
        "either could be considered correct, but they are the best options available."
    )
    print(textwrap.fill(conclusion_text, width=80))
    
    # In a multiple-choice setting where one answer must be selected, both C and D are excellent.
    # We will print all strong candidates identified by the script's logic.
    final_answer_text = f"The strongest candidates explaining both retention and divergence are choices {', '.join(best_choices)}."
    print(f"\n{final_answer_text}")


if __name__ == '__main__':
    analyze_gene_duplication_fates()