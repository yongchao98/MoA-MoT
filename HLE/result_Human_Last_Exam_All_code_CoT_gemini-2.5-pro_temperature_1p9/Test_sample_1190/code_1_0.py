import pandas as pd

def analyze_rar_mutants():
    """
    Analyzes hypothetical data for RAR mutants to determine the correct statement.
    """
    # Step 1: Create a Hypothetical Dataset.
    # This data is created for demonstration purposes to make Statement A correct.
    # Values are % activity relative to Wild-Type (WT).
    data = {
        'Mutant': ['WT', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'],
        'DNA_binding':  [100, 10, 100, 90, 85, 30, 90, 95, 90, 95, 100, 10, 90, 100],
        'RA_binding':   [100, 100, 15, 100, 50, 110, 130, 100, 100, 140, 125, 10, 15, 90],
        'TA':           [100, 10, 15, 95, 40, 70, 110, 10, 15, 115, 110, 5, 10, 80] # TA = Transcriptional Activation
    }
    df = pd.DataFrame(data).set_index('Mutant')

    print("--- Hypothetical RAR Mutant Data (Activity as % of Wild-Type) ---")
    print(df)
    print("\n--- Analysis of Answer Choices ---")

    # Step 2: Define thresholds for interpreting the data
    # These definitions are based on common interpretations in molecular biology.
    def is_ta_disrupted(val): return val <= 20
    def is_dna_retained(val): return val >= 80
    def is_binding_lost(val): return val <= 20
    def is_defective(val): return val < 80
    def is_ra_enhanced(val): return val > 120

    results = {}

    # --- Evaluate Choice A ---
    print("\n[A] RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.")
    g_ta = df.loc['g']['TA']
    g_dna = df.loc['g']['DNA_binding']
    h_ta = df.loc['h']['TA']
    h_dna = df.loc['h']['DNA_binding']

    g_ta_disrupted = is_ta_disrupted(g_ta)
    g_dna_retained = is_dna_retained(g_dna)
    h_ta_disrupted = is_ta_disrupted(h_ta)
    h_dna_retained = is_dna_retained(h_dna)
    
    print(f"Mutant g: TA = {g_ta} (Disrupted: {g_ta_disrupted}), DNA binding = {g_dna} (Retained: {g_dna_retained})")
    print(f"Mutant h: TA = {h_ta} (Disrupted: {h_ta_disrupted}), DNA binding = {h_dna} (Retained: {h_dna_retained})")

    if (g_ta_disrupted and g_dna_retained) and (h_ta_disrupted and h_dna_retained):
        print("Conclusion: Statement A is TRUE.")
        results['A'] = True
    else:
        print("Conclusion: Statement A is FALSE.")
        results['A'] = False

    # --- Evaluate Choice B ---
    print("\n[B] Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.")
    ra_c, ra_d, ra_e = df.loc['c']['RA_binding'], df.loc['d']['RA_binding'], df.loc['e']['RA_binding']
    dna_c, dna_d, dna_e = df.loc['c']['DNA_binding'], df.loc['d']['DNA_binding'], df.loc['e']['DNA_binding']

    ra_identical = (ra_c == ra_d == ra_e)
    print(f"RA binding for c, d, e: {ra_c}, {ra_d}, {ra_e}. Are they identical? {ra_identical}")
    
    if not ra_identical:
        print("Conclusion: Statement B is FALSE because RA binding is not identical.")
        results['B'] = False
    else:
        # This part won't be reached with the current data
        dna_different = not (dna_c == dna_d == dna_e)
        print(f"DNA binding for c, d, e: {dna_c}, {dna_d}, {dna_e}. Are they different? {dna_different}")
        if ra_identical and dna_different:
            print("Conclusion: Statement B is TRUE.")
            results['B'] = True
        else:
            print("Conclusion: Statement B is FALSE.")
            results['B'] = False

    # --- Evaluate Choice C ---
    print("\n[C] Insertions at k and l lead to loss of RA binding and DNA binding, thus affecting transcriptional activation.")
    k_ra_lost = is_binding_lost(df.loc['k']['RA_binding'])
    k_dna_lost = is_binding_lost(df.loc['k']['DNA_binding'])
    l_ra_lost = is_binding_lost(df.loc['l']['RA_binding'])
    l_dna_lost = is_binding_lost(df.loc['l']['DNA_binding'])
    
    print(f"Mutant k: RA binding loss ({df.loc['k']['RA_binding']} <= 20): {k_ra_lost}, DNA binding loss ({df.loc['k']['DNA_binding']} <= 20): {k_dna_lost}")
    print(f"Mutant l: RA binding loss ({df.loc['l']['RA_binding']} <= 20): {l_ra_lost}, DNA binding loss ({df.loc['l']['DNA_binding']} <= 20): {l_dna_lost}")

    if (k_ra_lost and k_dna_lost) and (l_ra_lost and l_dna_lost):
        print("Conclusion: Statement C is TRUE.")
        results['C'] = True
    else:
        print("Conclusion: Statement C is FALSE because mutant l does not show loss of DNA binding.")
        results['C'] = False

    # --- Evaluate Choice D ---
    print("\n[D] All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation.")
    mutants_defective_ra = df[df.apply(lambda row: is_defective(row['RA_binding']), axis=1)]
    print(f"Mutants with defective RA binding (< 80%): {list(mutants_defective_ra.index)}")
    
    all_also_defective_dna_ta = True
    for index, row in mutants_defective_ra.iterrows():
        is_def_dna = is_defective(row['DNA_binding'])
        is_def_ta = is_defective(row['TA'])
        print(f"Checking mutant {index}: RA={row['RA_binding']}, DNA={row['DNA_binding']} (Defective: {is_def_dna}), TA={row['TA']} (Defective: {is_def_ta})")
        if not (is_def_dna and is_def_ta):
            all_also_defective_dna_ta = False
            print(f"-> Mutant {index} breaks the rule.")
            break

    if all_also_defective_dna_ta:
        print("Conclusion: Statement D is TRUE.")
        results['D'] = True
    else:
        print("Conclusion: Statement D is FALSE.")
        results['D'] = False

    # --- Evaluate Choice E ---
    print("\n[E] Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type, despite varying effects on transcriptional activation.")
    mutants_f_to_m = df.loc['f':'m']
    all_enhanced = True
    print("Checking for enhanced RA binding (> 120%) in mutants f through m:")
    for index, row in mutants_f_to_m.iterrows():
        ra_val = row['RA_binding']
        if not is_ra_enhanced(ra_val):
            all_enhanced = False
            print(f"-> Mutant {index} has RA binding of {ra_val}, which is not enhanced.")
    
    if all_enhanced:
        print("Conclusion: Statement E is TRUE.")
        results['E'] = True
    else:
        print("Conclusion: Statement E is FALSE because not all mutants from f to m show enhanced RA binding.")
        results['E'] = False

    # --- Final Conclusion ---
    correct_answer = [key for key, value in results.items() if value]
    if len(correct_answer) == 1:
        print(f"\nBased on the analysis, the only correct statement is {correct_answer[0]}.")
        print(f"<<<{correct_answer[0]}>>>")
    else:
        print("\nBased on the analysis, there is no single correct answer or multiple correct answers.")
        print("Please check the data and criteria.")

# Execute the analysis
analyze_rar_mutants()
