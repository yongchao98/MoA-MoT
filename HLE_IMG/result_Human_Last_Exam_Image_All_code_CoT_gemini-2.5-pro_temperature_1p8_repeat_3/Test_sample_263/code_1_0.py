def analyze_emitter_stability():
    """
    Analyzes the stability of three Ir(III) complexes for use in LECs.
    """
    # Define the key features of each complex relevant to stability
    complexes = {
        'Complex 1': {
            'C^N Ligand': 'phenylpyridine (ppy)',
            'N^N Ligand': 'bipyridine (bpy)',
            'Features': ['Standard benchmark structure']
        },
        'Complex 2': {
            'C^N Ligand': 'phenylpyridine (ppy)',
            'N^N Ligand': 'large phenanthroimidazole derivative',
            'Features': ['Large, rigid ancillary ligand']
        },
        'Complex 3': {
            'C^N Ligand': '2,4-difluorophenylpyridine (dfppy)',
            'N^N Ligand': '4,4\'-di-tert-butyl-bipyridine (dtbbpy)',
            'Features': [
                'Fluorination on C^N ligand',
                'Bulky tert-butyl groups on N^N ligand'
            ]
        }
    }

    print("Analyzing the stability of Ir(III) emitters for LECs...\n")

    stability_assessment = {}
    for name, data in complexes.items():
        print(f"--- {name} ---")
        print(f"C^N Ligand: {data['C^N Ligand']}")
        print(f"N^N Ligand: {data['N^N Ligand']}")
        
        has_fluorination = 'fluoro' in data['C^N Ligand']
        has_bulk = 'tert-butyl' in data['N^N Ligand']
        
        reasons = []
        if has_fluorination:
            reasons.append("Fluorination increases resistance to oxidative degradation.")
        if has_bulk:
            reasons.append("Bulky groups provide steric protection against intermolecular degradation pathways.")
            
        if not reasons:
            reasons.append("This is a reference complex with no special stability-enhancing modifications.")
        
        stability_assessment[name] = reasons
        
        print("Stability Assessment:", " ".join(reasons))
        print("-" * (len(name) + 6) + "\n")

    print("--- Conclusion ---")
    print("Complex 3 includes two key features known to enhance device stability:")
    print("1. Fluorination of the main ligands, which improves chemical stability against oxidation.")
    print("2. Bulky tert-butyl groups on the ancillary ligand, which provide steric shielding.")
    print("\nTherefore, LECs based on Complex 3 are expected to be the most stable.")
    print("The correct answer choice is C.")

# Run the analysis
analyze_emitter_stability()