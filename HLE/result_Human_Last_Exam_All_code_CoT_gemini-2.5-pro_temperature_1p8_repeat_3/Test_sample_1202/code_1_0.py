import pandas as pd

def analyze_gene_duplication_fates():
    """
    Analyzes potential mechanisms for the retention and divergence of duplicate genes
    to determine the most likely candidate.
    """

    # Step 1: Define the mechanisms and their properties
    data = {
        'A': {
            'name': 'Gene conversion',
            'is_molecular_mechanism': True,
            'retains_both_functional': False,
            'causes_divergence': False,
            'explanation': 'Makes genes identical, preventing divergence.'
        },
        'B': {
            'name': 'Pseudogenization',
            'is_molecular_mechanism': True,
            'retains_both_functional': False,
            'causes_divergence': True, # Diverges into non-functionality
            'explanation': 'Leads to loss of one functional copy, not retention of two.'
        },
        'C': {
            'name': 'Neofunctionalization',
            'is_molecular_mechanism': True,
            'retains_both_functional': True,
            'causes_divergence': True,
            'explanation': 'One copy gains a new function; both are retained and diverge.'
        },
        'D': {
            'name': 'Subfunctionalization',
            'is_molecular_mechanism': True,
            'retains_both_functional': True,
            'causes_divergence': True,
            'explanation': 'Copies partition the original functions; both are required, retained, and diverge.'
        },
        'E': {
            'name': 'Adaptive radiation',
            'is_molecular_mechanism': False,
            'retains_both_functional': False,
            'causes_divergence': False,
            'explanation': 'A macro-evolutionary pattern, not a molecular mechanism at the gene level.'
        }
    }

    # For display purposes
    df = pd.DataFrame(data).T
    print("--- Analysis of Answer Choices ---")
    print(df[['name', 'is_molecular_mechanism', 'retains_both_functional', 'causes_divergence']])
    print("\n")

    # Step 2: Filter for mechanisms that meet the core requirements of the question.
    # The question asks for a mechanism that is responsible for BOTH retention and divergence.
    print("--- Filtering for Valid Mechanisms ---")
    print("Criterion 1: Must be a molecular mechanism acting on genes.")
    print("Criterion 2: Must lead to the retention of two functional gene copies.")
    print("Criterion 3: Must cause the gene copies to diverge in sequence/function.")

    valid_candidates = {}
    for choice, props in data.items():
        if props['is_molecular_mechanism'] and props['retains_both_functional'] and props['causes_divergence']:
            valid_candidates[choice] = props['name']

    print(f"\nMechanisms meeting all criteria: {valid_candidates}\n")

    # Step 3: Apply tie-breaking logic to determine the "most likely" mechanism.
    print("--- Tie-Breaker Analysis: Neofunctionalization vs. Subfunctionalization ---")
    print("Both Neofunctionalization (C) and Subfunctionalization (D) are valid and important models.")
    print("To determine the 'most likely', we consider the underlying requirements:")
    print(" - Neofunctionalization requires a relatively rare beneficial mutation for a new function to arise.")
    print(" - Subfunctionalization can occur via more common, neutral degenerative mutations, partitioning the ancestral functions.")
    print("\nBecause subfunctionalization provides a more passive and frequent pathway to preserve both copies from being lost, it is often considered a highly probable initial mechanism for retaining duplicate genes, which then allows them to diverge.")

    # Final conclusion
    final_answer_key = 'D' # Chosen based on the tie-breaker logic.
    final_answer_name = data[final_answer_key]['name']
    
    print("\n--- Conclusion ---")
    print(f"Based on this analysis, the most likely mechanism is ({final_answer_key}) {final_answer_name}.")

analyze_gene_duplication_fates()

print("\n<<<D>>>")