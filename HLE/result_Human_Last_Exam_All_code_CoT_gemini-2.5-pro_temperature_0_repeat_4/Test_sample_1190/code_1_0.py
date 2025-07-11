import pandas as pd
import io

def analyze_rar_mutants():
    """
    Analyzes hypothetical data for RAR mutants to determine the correct statement.
    Since the actual data was not provided, this function simulates the analysis
    with a plausible dataset.
    """
    # Step 1: Create a hypothetical dataset as a string in CSV format.
    # Values are percentages relative to Wild-Type (WT).
    hypothetical_data = """
    mutant,ra_binding,dna_binding,trans_activation
    WT,100,100,100
    a,90,95,85
    b,105,98,90
    c,50,80,40
    d,60,20,10
    e,50,90,45
    f,110,102,15
    g,95,105,10
    h,100,98,5
    i,99,40,20
    j,101,55,30
    k,5,10,2
    l,80,5,1
    m,115,108,25
    p,40,95,50
    """

    # Step 2: Load the data into a pandas DataFrame.
    df = pd.read_csv(io.StringIO(hypothetical_data)).set_index('mutant')
    print("--- Hypothetical Mutant Data Analysis ---")
    print(df)
    print("\n--- Evaluating Answer Choices ---")

    # Step 3: Define thresholds for analysis.
    disrupted_threshold = 20  # e.g., < 20% is disrupted/loss
    retained_threshold = 80   # e.g., > 80% is retained
    enhanced_threshold = 100  # > 100% is enhanced

    # --- Analysis for Choice A ---
    print("\n[Choice A]: RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.")
    mutant_g = df.loc['g']
    mutant_h = df.loc['h']
    g_trans_disrupted = mutant_g['trans_activation'] < disrupted_threshold
    h_trans_disrupted = mutant_h['trans_activation'] < disrupted_threshold
    g_dna_retained = mutant_g['dna_binding'] > retained_threshold
    h_dna_retained = mutant_h['dna_binding'] > retained_threshold
    
    is_A_true = g_trans_disrupted and h_trans_disrupted and g_dna_retained and h_dna_retained
    print(f"Mutant g: Trans. Activation = {mutant_g['trans_activation']} (Disrupted: {g_trans_disrupted}), DNA Binding = {mutant_g['dna_binding']} (Retained: {g_dna_retained})")
    print(f"Mutant h: Trans. Activation = {mutant_h['trans_activation']} (Disrupted: {h_trans_disrupted}), DNA Binding = {mutant_h['dna_binding']} (Retained: {h_dna_retained})")
    print(f"Conclusion for A: {is_A_true}")

    # --- Analysis for Choice B ---
    print("\n[Choice B]: Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.")
    mutants_b_group = df.loc[['c', 'd', 'e']]
    ra_binding_identical = mutants_b_group['ra_binding'].nunique() == 1
    dna_binding_differs = mutants_b_group['dna_binding'].std() > 20 # Check for significant variation
    
    is_B_true = ra_binding_identical and dna_binding_differs
    print(f"RA binding for c, d, e: {list(mutants_b_group['ra_binding'])}. Are they identical? {ra_binding_identical}")
    print(f"DNA binding for c, d, e: {list(mutants_b_group['dna_binding'])}. Do they differ significantly? {dna_binding_differs}")
    print(f"Conclusion for B: {is_B_true}")

    # --- Analysis for Choice C ---
    print("\n[Choice C]: Insertions at k and l lead to loss of RA binding and DNA binding, thus affecting transcriptional activation.")
    mutant_k = df.loc['k']
    mutant_l = df.loc['l']
    k_loss_ra = mutant_k['ra_binding'] < disrupted_threshold
    k_loss_dna = mutant_k['dna_binding'] < disrupted_threshold
    l_loss_ra = mutant_l['ra_binding'] < disrupted_threshold
    l_loss_dna = mutant_l['dna_binding'] < disrupted_threshold

    is_C_true = k_loss_ra and k_loss_dna and l_loss_ra and l_loss_dna
    print(f"Mutant k: RA binding loss? {k_loss_ra} (Value: {mutant_k['ra_binding']}), DNA binding loss? {k_loss_dna} (Value: {mutant_k['dna_binding']})")
    print(f"Mutant l: RA binding loss? {l_loss_ra} (Value: {mutant_l['ra_binding']}), DNA binding loss? {l_loss_dna} (Value: {mutant_l['dna_binding']})")
    print(f"Conclusion for C: {is_C_true}")

    # --- Analysis for Choice D ---
    print("\n[Choice D]: All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation.")
    defective_ra_mutants = df[df['ra_binding'] < 50]
    # Check if ALL of these are also defective in the other two categories
    all_also_defective = (defective_ra_mutants['dna_binding'] < 50).all() and (defective_ra_mutants['trans_activation'] < 50).all()
    
    is_D_true = all_also_defective
    print(f"Mutants defective in RA binding (<50%): {list(defective_ra_mutants.index)}")
    print("Checking if these are ALL also defective in DNA binding and Trans. Activation...")
    print(defective_ra_mutants[['dna_binding', 'trans_activation']])
    print(f"Conclusion for D: {is_D_true}")

    # --- Analysis for Choice E ---
    print("\n[Choice E]: Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type.")
    mutants_e_group = df.loc['f':'m']
    uniformly_enhanced = (mutants_e_group['ra_binding'] > enhanced_threshold).all()

    is_E_true = uniformly_enhanced
    print(f"RA binding for mutants f through m: {list(mutants_e_group['ra_binding'])}")
    print(f"Are they all enhanced (>100)? {uniformly_enhanced}")
    print(f"Conclusion for E: {is_E_true}")

if __name__ == '__main__':
    analyze_rar_mutants()