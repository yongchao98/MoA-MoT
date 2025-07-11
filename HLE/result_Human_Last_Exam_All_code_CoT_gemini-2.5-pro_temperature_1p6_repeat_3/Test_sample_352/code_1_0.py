import pandas as pd
from io import StringIO

def suggest_mutagenesis():
    """
    Suggests and explains the site-directed mutagenesis plan for protein x
    to relieve the inhibitory, negatively charged patch.
    """
    # Define the original amino acids and their properties
    original_aa_data = {
        'Position': [47, 48, 49, 50],
        'AA_Code': ['S', 'E', 'E', 'D'],
        'Name': ['Serine', 'Glutamate', 'Glutamate', 'Aspartate'],
        'Property': [
            'Polar, phosphorylatable (-OH group)',
            'Acidic, negatively charged',
            'Acidic, negatively charged',
            'Acidic, negatively charged'
        ]
    }
    
    # Define the replacement amino acid and its properties
    replacement_aa = {
        'code': 'A',
        'name': 'Alanine',
        'property': 'Non-polar, neutral charge, small size',
        'rationale': (
            "This replacement neutralizes the patch by removing negative charges "
            "and prevents phosphorylation at position 47. Its small size minimizes "
            "potential structural disruption, making it the ideal choice."
        )
    }

    print("--- Site-Directed Mutagenesis Plan ---")
    print("\nObjective: To eliminate the negative charge from the S-E-E-D patch (positions 47-50).")
    print(f"Recommendation: Replace all four residues with {replacement_aa['name']} ({replacement_aa['code']}).")
    print(f"\nRationale: {replacement_aa['rationale']}\n")
    print("--- Mutation Details ---")

    original_sequence = []
    final_sequence = []
    
    # Print the detailed plan for each position
    for i, pos in enumerate(original_aa_data['Position']):
        original_code = original_aa_data['AA_Code'][i]
        original_name = original_aa_data['Name'][i]
        original_prop = original_aa_data['Property'][i]
        
        print(f"Position {pos}:")
        print(f"  - Original:   {original_name} ({original_code}) | Property: {original_prop}")
        print(f"  - Proposed:   {replacement_aa['name']} ({replacement_aa['code']}) | Property: {replacement_aa['property']}")
        print(f"  - Mutation:   {original_code}{pos}{replacement_aa['code']}")
        
        original_sequence.append(original_code)
        final_sequence.append(replacement_aa['code'])
    
    # Print the final summary equation
    original_str = "".join(original_sequence)
    final_str = "".join(final_sequence)
    
    print("\n--- Final Equation of Change ---")
    print(f"Original sequence at 47-50: {original_str}")
    print(f"Mutated sequence at 47-50:  {final_str}")
    print(f"\nOverall Change: {original_str} -> {final_str}")


suggest_mutagenesis()
<<<AAAA>>>