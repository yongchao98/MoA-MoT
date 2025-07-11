import pandas as pd

def design_mutagenesis_experiment():
    """
    Designs and explains the site-directed mutagenesis experiment to neutralize
    the negatively charged patch in protein x.
    """
    # Define the original and replacement amino acids
    original_patch = {
        47: {'name': 'Serine', 'code': 'S', 'charge': 'Neutral (becomes negative upon phosphorylation)'},
        48: {'name': 'Glutamate', 'code': 'E', 'charge': 'Negative'},
        49: {'name': 'Glutamate', 'code': 'E', 'charge': 'Negative'},
        50: {'name': 'Aspartate', 'code': 'D', 'charge': 'Negative'}
    }
    
    replacement_amino_acid = {'name': 'Alanine', 'code': 'A', 'charge': 'Neutral'}

    print("### Site-Directed Mutagenesis Plan for Protein X ###\n")
    print("Objective: To eliminate the autoinhibitory negative charge from the IDR patch (residues 47-50).\n")
    
    print("Rationale:")
    print("The S-E-E-D patch is hypothesized to be inhibitory due to its strong negative charge. To test this, we will replace all four residues with Alanine (A). Alanine is small, neutral, and cannot be phosphorylated, making it the ideal choice to neutralize the patch without introducing other significant chemical or structural changes.\n")

    print("Proposed Mutations:")
    
    mutations_list = []
    for position, aa_info in original_patch.items():
        original_aa = f"{aa_info['name']} ({aa_info['code']})"
        replacement_aa = f"{replacement_amino_acid['name']} ({replacement_amino_acid['code']})"
        mutation_str = f"Position {position}: {original_aa} -> {replacement_aa}"
        print(mutation_str)
        mutations_list.append(f"{aa_info['code']}{position}{replacement_amino_acid['code']}")

    final_mutant_sequence = "".join([replacement_amino_acid['code']] * len(original_patch))

    print("\n-----------------------------------------------------")
    print(f"Original Sequence (47-50): {'-'.join([v['code'] for v in original_patch.values()])}")
    print(f"Proposed Mutant Sequence (47-50): {'-'.join(list(final_mutant_sequence))}")
    print(f"Standard Mutation Notation: {', '.join(mutations_list)}")
    print("-----------------------------------------------------")

if __name__ == '__main__':
    design_mutagenesis_experiment()
