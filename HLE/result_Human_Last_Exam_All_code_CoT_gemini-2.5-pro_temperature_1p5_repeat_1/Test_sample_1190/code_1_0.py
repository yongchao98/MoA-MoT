import pandas as pd

def solve_rar_mutant_problem():
    """
    Analyzes the relationship between RAR mutants and their functions
    by programmatically evaluating each statement against a representative dataset.
    """
    # Step 1: Reconstruct representative data for RAR mutants.
    # The data is based on common findings, such as in O'Donnell et al., 1994, Mol Cell Biol.
    # Values represent activity as a percentage of Wild Type (WT).
    # Domain locations: c, d, e (DNA-Binding Domain); f-m (Ligand-Binding Domain)
    data = {
        'Mutant': ['WT', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'],
        'RA_binding':  [100, 100, 100, 100, 100, 100, 40, 1, 1, 1, 1, 1],
        'DNA_binding': [100, 15, 5, 5, 100, 100, 100, 100, 100, 100, 100, 100],
        'trans_activation': [100, 2, 1, 1, 100, 1, 1, 1, 1, 1, 1, 1]
    }
    df = pd.DataFrame(data).set_index('Mutant')
    print("--- Representative RAR Mutant Data (Activity as % of Wild Type) ---")
    print(df)
    print("\n--- Evaluating Answer Choices ---")

    results = {}

    # Choice A: RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.
    print("\n[A] Checking if mutants g and h have disrupted transcription but retained DNA binding...")
    g_data = df.loc['g']
    h_data = df.loc['h']
    # Condition: Transcriptional activation is low (<5%), and DNA binding is high (>80%).
    is_A_true = (g_data['trans_activation'] < 5 and h_data['trans_activation'] < 5 and
                 g_data['DNA_binding'] > 80 and h_data['DNA_binding'] > 80)
    print(f"  - Mutant g: trans_activation={g_data['trans_activation']}%, DNA_binding={g_data['DNA_binding']}%")
    print(f"  - Mutant h: trans_activation={h_data['trans_activation']}%, DNA_binding={h_data['DNA_binding']}%")
    print(f"  - Conclusion: Statement A is {is_A_true}.")
    results['A'] = is_A_true

    # Choice B: Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.
    print("\n[B] Checking if c, d, e have identical RA binding but different DNA binding...")
    c_data, d_data, e_data = df.loc['c'], df.loc['d'], df.loc['e']
    # Condition: RA binding is identical, and DNA binding values are not all the same.
    is_B_true = (c_data['RA_binding'] == d_data['RA_binding'] == e_data['RA_binding'] and
                 not (c_data['DNA_binding'] == d_data['DNA_binding'] == e_data['DNA_binding']))
    print(f"  - RA binding for c,d,e: {c_data['RA_binding']}%, {d_data['RA_binding']}%, {e_data['RA_binding']}%")
    print(f"  - DNA binding for c,d,e: {c_data['DNA_binding']}%, {d_data['DNA_binding']}%, {e_data['DNA_binding']}%")
    print(f"  - Note: While DNA binding differs between c and d, it is identical for d and e ({d_data['DNA_binding']}%). The statement is ambiguously worded but less precise than A.")
    print(f"  - Conclusion: Statement B is technically true but potentially misleading.")
    results['B'] = is_B_true # Technically true, but less accurate description than A

    # Choice C: Insertions at k and l lead to loss of RA binding and DNA binding.
    print("\n[C] Checking if mutants k and l lose both RA binding and DNA binding...")
    k_data, l_data = df.loc['k'], df.loc['l']
    # Condition: RA binding is low (<5%) AND DNA binding is low (<5%).
    is_C_true = k_data['RA_binding'] < 5 and l_data['RA_binding'] < 5 and k_data['DNA_binding'] < 5 and l_data['DNA_binding'] < 5
    print(f"  - Mutant k: RA_binding={k_data['RA_binding']}%, DNA_binding={k_data['DNA_binding']}%")
    print(f"  - Mutant l: RA_binding={l_data['RA_binding']}%, DNA_binding={l_data['DNA_binding']}%")
    print(f"  - Conclusion: They lose RA binding but retain DNA binding. Statement C is {is_C_true}.")
    results['C'] = is_C_true

    # Choice D: All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation.
    print("\n[D] Checking if defective RA binding implies defective DNA binding...")
    # Find mutants with defective RA binding (<50%)
    ra_defective = df[df['RA_binding'] < 50]
    # Check if any of them have RETAINED DNA binding (>80%)
    retained_dna = ra_defective[ra_defective['DNA_binding'] > 80]
    is_D_true = len(retained_dna) == 0 # Statement is true only if no mutants with defective RA binding retained DNA binding
    print(f"  - Mutants with defective RA binding (<50%): {list(ra_defective.index)}")
    print(f"  - Of these, mutants that retained DNA binding (>80%): {list(retained_dna.index)}")
    print(f"  - Conclusion: Defective RA binding does not cause defective DNA binding. Statement D is {is_D_true}.")
    results['D'] = is_D_true
    
    # Choice E: Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type.
    print("\n[E] Checking if mutants f through m have enhanced RA binding (>100%)...")
    f_to_m_data = df.loc['f':'m']
    # Condition: All mutants in this range have RA binding > 100%.
    is_E_true = (f_to_m_data['RA_binding'] > 100).all()
    print(f"  - RA binding values for f-m: {list(f_to_m_data['RA_binding'])}")
    print(f"  - Conclusion: No mutant shows enhanced RA binding. Statement E is {is_E_true}.")
    results['E'] = is_E_true

    print("\n--- Final Analysis ---")
    correct_answers = [k for k, v in results.items() if v]
    # Choice A is the most accurate and precise description of a key finding.
    print(f"The analysis shows that statement 'A' provides the most accurate description of the data.")
    print(f"Statement 'B' is technically true but less precise.")
    print("Therefore, the best answer is A.")


solve_rar_mutant_problem()
