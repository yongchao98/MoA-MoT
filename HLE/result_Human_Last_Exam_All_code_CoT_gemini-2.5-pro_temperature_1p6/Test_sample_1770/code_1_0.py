def find_knockout_target():
    """
    Identifies the optimal gene for knockout to prevent p-coumaric acid degradation
    by searching a simplified mock database for Corynebacterium glutamicum.
    """
    # A simplified mock database of C. glutamicum genes related to p-coumaric acid.
    # In a real-world scenario, this data would come from databases like KEGG, BioCyc, or NCBI.
    c_glutamicum_db = {
        'tal': {
            'enzyme': 'Tyrosine ammonia-lyase',
            'reaction': 'L-tyrosine -> p-coumaric acid + NH3',
            'function': 'Biosynthesis of p-coumaric acid'
        },
        'fcs': {
            'locus_tag': 'cg1180',
            'enzyme': 'Feruloyl-CoA synthetase / p-Coumarate-CoA ligase',
            'reaction': 'p-coumaric acid + CoA + ATP -> p-coumaroyl-CoA + AMP + PPi',
            'function': 'Degradation of p-coumaric acid (first step)'
        },
        'ech': {
            'locus_tag': 'cg1179',
            'enzyme': 'Enoyl-CoA hydratase / P-coumaroyl-CoA hydratase/lyase',
            'reaction': 'p-coumaroyl-CoA -> 4-hydroxybenzoyl-CoA + Acetyl-CoA',
            'function': 'Degradation of p-coumaric acid (second step)'
        },
        'pyk': {
            'enzyme': 'Pyruvate kinase',
            'reaction': 'Phosphoenolpyruvate -> Pyruvate',
            'function': 'Central metabolism (Glycolysis)'
        }
    }
    
    target_compound = 'p-coumaric acid'
    process_to_block = 'Degradation'
    
    print(f"Searching for gene to knock out to prevent {process_to_block} of {target_compound}...\n")
    
    best_target = None

    # Search the database for the gene that catalyzes the first step of degradation
    for gene, data in c_glutamicum_db.items():
        is_degradation_gene = process_to_block.lower() in data.get('function', '').lower()
        is_first_step = 'first step' in data.get('function', '').lower()
        
        # Split the reaction to identify substrates
        if 'reaction' in data:
            substrates = data['reaction'].split(' -> ')[0]
            is_correct_substrate = target_compound in substrates
        else:
            is_correct_substrate = False

        if is_degradation_gene and is_first_step and is_correct_substrate:
            best_target = (gene, data)
            break
            
    if best_target:
        gene_name, gene_data = best_target
        print("--- Knockout Target Identified ---")
        print(f"Gene to knock out: {gene_name}")
        if 'locus_tag' in gene_data:
            print(f"Locus Tag (in C. glutamicum ATCC 13032): {gene_data['locus_tag']}")
        print(f"Enzyme: {gene_data['enzyme']}")
        print(f"Reaction Catalyzed: {gene_data['reaction']}")
        print("\nReasoning:")
        print(f"The '{gene_name}' gene encodes the enzyme that catalyzes the first committed step in the degradation of {target_compound}.")
        print("By knocking out this gene, p-coumaric acid can no longer enter the degradation pathway and will accumulate in the cell.")
    else:
        print("Could not identify an optimal knockout target from the database.")

if __name__ == '__main__':
    find_knockout_target()
