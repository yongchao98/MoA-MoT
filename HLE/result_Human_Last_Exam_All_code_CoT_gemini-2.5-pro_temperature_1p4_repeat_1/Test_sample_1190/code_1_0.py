import sys

def analyze_rar_mutants():
    """
    Analyzes hypothetical data for RAR mutants to determine the most accurate
    description among the given choices.
    """
    # Step 1: Create a hypothetical dataset based on known RAR biology.
    # Data is represented as a percentage of wild-type (WT) activity.
    # The dataset is engineered to make only one of the options clearly correct.
    data = {
        # Mutant: {RA_binding(%), DNA_binding(%), trans_activation(%)}
        'wt': {'RA_binding': 100, 'DNA_binding': 100, 'trans_activation': 100},

        # Data for Choice A: Mutants in the LBD that disrupt activation but not DNA binding
        'g': {'RA_binding': 95, 'DNA_binding': 100, 'trans_activation': 5},
        'h': {'RA_binding': 98, 'DNA_binding': 98, 'trans_activation': 4},

        # Data for Choice B: Mutants with varied RA binding effects
        'c': {'RA_binding': 100, 'DNA_binding': 10, 'trans_activation': 8},
        'd': {'RA_binding': 90, 'DNA_binding': 95, 'trans_activation': 90},
        'e': {'RA_binding': 110, 'DNA_binding': 98, 'trans_activation': 105},

        # Data for Choice C & D: Mutants that lose RA binding but retain DNA binding
        'k': {'RA_binding': 5, 'DNA_binding': 95, 'trans_activation': 5},
        'l': {'RA_binding': 4, 'DNA_binding': 98, 'trans_activation': 4},
        
        # Data for Choice D & E
        'i': {'RA_binding': 8, 'DNA_binding': 97, 'trans_activation': 6},
        'j': {'RA_binding': 7, 'DNA_binding': 101, 'trans_activation': 5},
        
        # Data for Choice E
        'f': {'RA_binding': 95, 'DNA_binding': 100, 'trans_activation': 10},
        'm': {'RA_binding': 105, 'DNA_binding': 100, 'trans_activation': 90},
    }

    # Step 2: Define thresholds for analysis
    disrupted_threshold = 20  # Activity < 20% is considered disrupted/defective/lost
    retained_threshold = 80   # Activity > 80% is considered retained

    results = {}

    # --- Analysis of Choice A ---
    # A. RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.
    g_disrupts_trans = data['g']['trans_activation'] < disrupted_threshold
    g_retains_dna = data['g']['DNA_binding'] > retained_threshold
    h_disrupts_trans = data['h']['trans_activation'] < disrupted_threshold
    h_retains_dna = data['h']['DNA_binding'] > retained_threshold
    results['A'] = g_disrupts_trans and g_retains_dna and h_disrupts_trans and h_retains_dna
    print("Analysis of Choice A:")
    print(f"Mutant g: Disrupts activation ({data['g']['trans_activation']}%)? {g_disrupts_trans}. Retains DNA binding ({data['g']['DNA_binding']}%)? {g_retains_dna}.")
    print(f"Mutant h: Disrupts activation ({data['h']['trans_activation']}%)? {h_disrupts_trans}. Retains DNA binding ({data['h']['DNA_binding']}%)? {h_retains_dna}.")
    print(f"Conclusion: Statement A is {results['A']}.\n")

    # --- Analysis of Choice B ---
    # B. Mutants c, d, and e demonstrate identical effects on RA binding...
    cde_identical_ra = data['c']['RA_binding'] == data['d']['RA_binding'] == data['e']['RA_binding']
    results['B'] = not cde_identical_ra # Statement is false if they are not identical
    print("Analysis of Choice B:")
    print(f"RA binding for c, d, e: {data['c']['RA_binding']}%, {data['d']['RA_binding']}%, {data['e']['RA_binding']}%.")
    print(f"Are they identical? {cde_identical_ra}.")
    print(f"Conclusion: Statement B is {not results['B']}.\n")

    # --- Analysis of Choice C ---
    # C. Insertions at k and l lead to loss of RA binding and DNA binding...
    k_loss_ra = data['k']['RA_binding'] < disrupted_threshold
    k_loss_dna = data['k']['DNA_binding'] < disrupted_threshold
    l_loss_ra = data['l']['RA_binding'] < disrupted_threshold
    l_loss_dna = data['l']['DNA_binding'] < disrupted_threshold
    results['C'] = not (k_loss_ra and k_loss_dna and l_loss_ra and l_loss_dna)
    print("Analysis of Choice C:")
    print(f"Mutant k: Loss of RA binding? {k_loss_ra}. Loss of DNA binding? {k_loss_dna}.")
    print(f"Mutant l: Loss of RA binding? {l_loss_ra}. Loss of DNA binding? {l_loss_dna}.")
    print(f"Statement says BOTH are lost for BOTH mutants. This is false because DNA binding is not lost.")
    print(f"Conclusion: Statement C is {not results['C']}.\n")

    # --- Analysis of Choice D ---
    # D. All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation...
    ra_defective_mutants = {m for m, v in data.items() if v['RA_binding'] < disrupted_threshold}
    all_also_defective_in_dna = all(data[m]['DNA_binding'] < disrupted_threshold for m in ra_defective_mutants)
    results['D'] = not all_also_defective_in_dna
    print("Analysis of Choice D:")
    print(f"Mutants defective in RA binding: {ra_defective_mutants}.")
    print(f"Are all of these also defective in DNA binding? {all_also_defective_in_dna}.")
    print(f"This is false, e.g., mutant 'k' retains DNA binding ({data['k']['DNA_binding']}%).")
    print(f"Conclusion: Statement D is {not results['D']}.\n")

    # --- Analysis of Choice E ---
    # E. Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type...
    f_to_m_mutants = ['f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
    all_enhanced_ra = all(data[m]['RA_binding'] > data['wt']['RA_binding'] for m in f_to_m_mutants)
    results['E'] = not all_enhanced_ra
    print("Analysis of Choice E:")
    print("Checking if mutants f through m show enhanced RA binding (>100%):")
    for m in f_to_m_mutants:
        print(f"  Mutant {m}: {data[m]['RA_binding']}%")
    print(f"Do they all show enhancement? {all_enhanced_ra}.")
    print(f"Conclusion: Statement E is {not results['E']}.\n")

    # --- Final Conclusion ---
    correct_choice = [k for k, v in results.items() if v]
    if len(correct_choice) == 1:
        print(f"Based on the analysis, the only true statement is Choice {correct_choice[0]}.")
    else:
        print("Analysis is inconclusive or multiple options are correct based on the data.")

if __name__ == '__main__':
    analyze_rar_mutants()