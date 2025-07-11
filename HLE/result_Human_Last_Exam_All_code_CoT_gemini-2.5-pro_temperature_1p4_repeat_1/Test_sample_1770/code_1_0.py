import textwrap

def find_gene_for_knockout():
    """
    This script queries a simplified metabolic database for Corynebacterium glutamicum
    to identify the optimal gene for knockout to prevent product degradation.
    The specific goal is to stop the degradation of p-coumaric acid.
    """
    # A simplified database representing known degradation pathways.
    # In a real-world scenario, this data would be sourced from genomic and metabolic
    # databases like KEGG, BioCyc, or NCBI.
    c_glutamicum_degradation_db = {
        'p-coumaric acid': {
            'gene_name': 'fcs',
            'locus_tag': 'cg2904',
            'enzyme': 'p-coumaroyl-CoA synthetase (also known as Feruloyl-CoA synthetase)',
            'reaction': 'p-coumaric acid + CoA -> p-coumaroyl-CoA',
            'rationale': ('This enzyme catalyzes the very first step of the p-coumaric acid '
                          'degradation pathway. By attaching Coenzyme A (CoA) to p-coumaric acid, '
                          'it "activates" the molecule for subsequent breakdown. Deleting this '
                          'gene is the most direct and effective strategy to prevent the loss '
                          'of p-coumaric acid, thereby increasing the final product yield.')
        },
        'protocatechuate': {
            'gene_name': 'pcaGH',
            'enzyme': 'Protocatechuate 3,4-dioxygenase',
            'reaction': 'Protocatechuate -> Ring Cleavage',
            'rationale': ('This enzyme breaks the aromatic ring of protocatechuate, which is a '
                          'central intermediate in the degradation of many aromatic compounds. '
                          'While a knockout here would also stop degradation, it would lead to '
                          'the accumulation of the intermediate protocatechuate, not the desired '
                          'p-coumaric acid.')
        }
    }

    # The compound we want to produce and prevent from being degraded.
    target_compound = 'p-coumaric acid'
    
    # Query the database for the target compound.
    knockout_info = c_glutamicum_degradation_db.get(target_compound)

    print(f"Query: Identify gene to knock out to prevent degradation of '{target_compound}'.")
    print("-" * 75)

    if knockout_info:
        print(">>> Optimal Knockout Target Found <<<\n")
        print(f"Gene to Knock Out:         {knockout_info['gene_name']}")
        print(f"Encoded Enzyme:            {knockout_info['enzyme']}")
        print(f"Blocked Reaction:          {knockout_info['reaction']}")
        print("\nRationale:")
        
        # Use textwrap for cleaner formatting of the rationale text.
        wrapped_rationale = textwrap.fill(knockout_info['rationale'], width=75)
        print(wrapped_rationale)
        print("-" * 75)
    else:
        print(f"Could not find a recommended knockout target for {target_compound}.")

# Execute the function to find and print the answer.
find_gene_for_knockout()