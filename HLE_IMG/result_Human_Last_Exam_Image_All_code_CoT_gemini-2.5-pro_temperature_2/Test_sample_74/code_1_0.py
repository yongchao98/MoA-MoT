import numpy as np

def evaluate_statements():
    """
    Evaluates claims about the data in Figure 3b of Siegel et al. 2024.
    NGSC values are visually estimated from the graph for each participant.
    """
    
    # Visually estimated NGSC data from Figure 3b (right)
    # The data represents the approximate range and distribution of points.
    participants_data = {
        'P1': {
            'no_drug': [0.66, 0.68, 0.70, 0.71, 0.72],
            'mtp': [0.68, 0.69, 0.70, 0.71],
            'psilocybin': [0.71, 0.72, 0.73, 0.75]
        },
        'P3': {
            'no_drug': [0.63, 0.65, 0.68, 0.70, 0.73],
            'mtp': [0.66, 0.68, 0.69, 0.70],
            'psilocybin': [0.68, 0.69, 0.70, 0.72]
        },
        'P4': {
            'no_drug': [0.65, 0.67, 0.69, 0.71, 0.72],
            'mtp': [0.66, 0.68, 0.70],
            'psilocybin': [0.66, 0.75, 0.78, 0.82]
        },
        'P5': {
            'no_drug': [0.61, 0.63, 0.65, 0.68, 0.70],
            'mtp': [0.62, 0.64, 0.66],
            'psilocybin': [0.71, 0.73, 0.74, 0.76]
        },
        'P6': { # Data for P6 can also be seen clearly on the left panel
            'no_drug': [0.63, 0.65, 0.66, 0.67, 0.68, 0.69],
            'mtp': [0.64, 0.65, 0.66],
            'psilocybin': [0.72, 0.74]
        },
        'P7': {
            'no_drug': [0.65, 0.67, 0.69, 0.71, 0.72],
            'mtp': [0.66, 0.68, 0.70],
            'psilocybin': [0.76, 0.78, 0.79, 0.81]
        }
    }

    # --- Analysis of choice F ---
    # F. At least one participant had the lowest global NGSC measurement 
    #    during a scan under the psilocybin condition.
    
    print("--- Evaluating Statement F ---")
    is_statement_F_true = False
    
    for p_id, data in participants_data.items():
        all_scans = data['no_drug'] + data['mtp'] + data['psilocybin']
        min_ngsc_overall = min(all_scans)
        
        # Check if this minimum value belongs to the psilocybin scans
        if min_ngsc_overall in data['psilocybin']:
            is_statement_F_true = True
            print(f"For participant {p_id}, the lowest NGSC ({min_ngsc_overall}) was in the psilocybin condition. Statement F might be true.")
            break
        else:
            min_ngsc_psilo = min(data['psilocybin'])
            print(f"For participant {p_id}:")
            print(f"  Lowest NGSC overall = {min_ngsc_overall}")
            print(f"  Lowest NGSC in Psilocybin condition = {min_ngsc_psilo}")
            print("  The overall lowest NGSC was not from a psilocybin scan.")

    if not is_statement_F_true:
        print("\nConclusion for F: Statement F is false. For every participant shown, the minimum NGSC value occurs in a non-psilocybin condition.")

    print("\n--- Summary of Other Options ---")
    print("Option B: Unsupported. Requires inferring statistical significance which is not provided.")
    print("Option D: False. High entropy (even distribution) means low, not high, functional connectivity.")
    print("Option K: False. For P4, the lowest psilocybin NGSC (~0.66) is not higher than all other scans (e.g., a no-drug scan at ~0.72).")
    print("Other options are incorrect based on definitions (NGSC must be >= 0) or experimental design (not pre- vs post-).")
    print("\nSince options B-N are incorrect or unsupported, the correct option is A.")

evaluate_statements()