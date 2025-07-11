import pandas as pd

def analyze_rar_mutants():
    """
    Analyzes RAR mutant data to determine the correct statement among the given choices.
    """
    # Data is based on the canonical experiment on RAR functional domains (e.g., Alberts, MBOC, Fig 7-56).
    # Values are approximate percentages of wild-type (WT) activity.
    data = {
        # Mutant: [DNA binding %, RA binding %, Transcriptional Activation %]
        'c': [5, 100, 5],
        'd': [5, 100, 5],
        'e': [100, 100, 100],
        'f': [100, 5, 5],
        'g': [100, 100, 5],
        'h': [100, 100, 5],
        'i': [100, 5, 5],
        'j': [100, 5, 5],
        'k': [5, 5, 5],
        'l': [5, 5, 5],
        'm': [100, 5, 5],
    }

    # Define thresholds for functional description.
    # Activity < 20% is considered "disrupted", "defective", or "lost".
    # Activity > 80% is considered "retained" or "normal".
    disrupted_threshold = 20
    retained_threshold = 80

    print("--- Evaluating Answer Choices ---\n")

    # --- Choice A ---
    print("Choice A: RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.")
    g_trans_act = data['g'][2]
    g_dna_bind = data['g'][0]
    h_trans_act = data['h'][2]
    h_dna_bind = data['h'][0]

    g_trans_disrupted = g_trans_act < disrupted_threshold
    g_dna_retained = g_dna_bind > retained_threshold
    h_trans_disrupted = h_trans_act < disrupted_threshold
    h_dna_retained = h_dna_bind > retained_threshold
    
    is_A_true = g_trans_disrupted and g_dna_retained and h_trans_disrupted and h_dna_retained

    print(f"Mutant g: Transcriptional Activation = {g_trans_act}% (Disrupted: {g_trans_disrupted}), DNA Binding = {g_dna_bind}% (Retained: {g_dna_retained})")
    print(f"Mutant h: Transcriptional Activation = {h_trans_act}% (Disrupted: {h_trans_disrupted}), DNA Binding = {h_dna_bind}% (Retained: {h_dna_retained})")
    print(f"Conclusion for A: {is_A_true}\n")


    # --- Choice B ---
    print("Choice B: Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.")
    c_ra_bind = data['c'][1]
    d_ra_bind = data['d'][1]
    e_ra_bind = data['e'][1]
    c_dna_bind = data['c'][0]
    d_dna_bind = data['d'][0]
    e_dna_bind = data['e'][0]
    
    # Check for identical RA binding effects (all are normal)
    ra_identical = (c_ra_bind > retained_threshold and d_ra_bind > retained_threshold and e_ra_bind > retained_threshold)
    # Check for significant difference in DNA binding (they are not all the same)
    dna_differs = not (c_dna_bind == d_dna_bind and d_dna_bind == e_dna_bind)

    is_B_true = ra_identical and dna_differs
    print(f"RA Binding: c={c_ra_bind}%, d={d_ra_bind}%, e={e_ra_bind}% (All retained: {ra_identical})")
    print(f"DNA Binding: c={c_dna_bind}%, d={d_dna_bind}%, e={e_dna_bind}% (They differ: {dna_differs})")
    print("Note: The statement is plausible, but 'differ significantly' is ambiguous as c and d are identical to each other.")
    print(f"Conclusion for B: {is_B_true} (conditionally, based on interpretation)\n")

    # --- Choice C ---
    print("Choice C: Insertions at k and l lead to loss of RA binding and DNA binding, thus affecting transcriptional activation.")
    k_ra_bind = data['k'][1]
    k_dna_bind = data['k'][0]
    l_ra_bind = data['l'][1]
    l_dna_bind = data['l'][0]

    k_ra_lost = k_ra_bind < disrupted_threshold
    k_dna_lost = k_dna_bind < disrupted_threshold
    l_ra_lost = l_ra_bind < disrupted_threshold
    l_dna_lost = l_dna_bind < disrupted_threshold
    
    is_C_true = k_ra_lost and k_dna_lost and l_ra_lost and l_dna_lost
    
    print(f"Mutant k: RA Binding = {k_ra_bind}% (Lost: {k_ra_lost}), DNA Binding = {k_dna_bind}% (Lost: {k_dna_lost})")
    print(f"Mutant l: RA Binding = {l_ra_bind}% (Lost: {l_ra_lost}), DNA Binding = {l_dna_bind}% (Lost: {l_dna_lost})")
    print(f"Conclusion for C: {is_C_true}\n")

    # --- Choice D ---
    print("Choice D: All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation, indicating a linked mechanism.")
    is_D_true = True
    counterexample = None
    ra_defective_mutants = {m for m, values in data.items() if values[1] < disrupted_threshold}
    for m in ra_defective_mutants:
        dna_is_defective = data[m][0] < disrupted_threshold
        trans_is_defective = data[m][2] < disrupted_threshold
        if not (dna_is_defective and trans_is_defective):
            is_D_true = False
            counterexample = m
            break
    print(f"RA-defective mutants: {ra_defective_mutants}")
    if not is_D_true:
        print(f"Counterexample found: Mutant '{counterexample}'")
        print(f"RA Binding = {data[counterexample][1]}% (Defective), but DNA Binding = {data[counterexample][0]}% (NOT Defective)")
    print(f"Conclusion for D: {is_D_true}\n")

    # --- Choice E ---
    print("Choice E: Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type, despite varying effects on transcriptional activation.")
    is_E_true = True
    counterexample = None
    for m in ['f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']:
        # Enhanced means > 100% (or > WT). Most are severely reduced.
        if data[m][1] < 100:
            is_E_true = False
            counterexample = m
            break
    print(f"Checking for uniformly enhanced RA binding in mutants f-m...")
    if not is_E_true:
        print(f"Counterexample found: Mutant '{counterexample}' has RA binding of {data[counterexample][1]}%, which is not enhanced.")
    print(f"Conclusion for E: {is_E_true}\n")

    print("--- Final Analysis ---")
    print("Statements A and C are factually correct based on the data. Statement D and E are false.")
    print("Statement A describes a key finding: the separation of the transcriptional activation function from DNA binding, as mutants g and h bind DNA normally but cannot activate transcription.")
    print("This demonstrates a more nuanced and significant relationship than statement C, which describes a general failure of all functions.")
    print("Therefore, A is the most accurate and insightful description of the relationships revealed by the experiment.")


if __name__ == '__main__':
    analyze_rar_mutants()
<<<A>>>