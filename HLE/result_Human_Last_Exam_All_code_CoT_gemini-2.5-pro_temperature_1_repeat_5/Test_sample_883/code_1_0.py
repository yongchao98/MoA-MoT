def analyze_enzyme_data():
    """
    Analyzes experimental data for enzyme Zma1 to determine the function of
    various molecules and select the correct conclusion from a list of choices.
    """
    # Store the experimental results in a dictionary
    results = {
        'Control': 500,
        'MgCl2': 700,
        'CaCl2': 500,
        'CuCl2': 400,
        'Al1': 1000,
        'Al2': 150,
        'Al1_Al2': 150,
        'XAG1': 10,
        'XAG1_plus_A': 450,
        'Rga1': 10,
        'Rga1_plus_A': 10,
    }

    kcat_control = results['Control']

    print("--- Analysis of Molecule Al1 ---")
    kcat_al1 = results['Al1']
    print(f"The kcat for the control reaction is {kcat_control}/second.")
    print(f"When Al1 is added, the kcat changes to {kcat_al1}/second.")
    if kcat_al1 > kcat_control:
        print("Conclusion: Since the rate increased, Al1 functions as an activator.\n")
    else:
        print("Conclusion: Al1 is not an activator.\n")

    print("--- Analysis of Molecule Rga1 ---")
    kcat_rga1 = results['Rga1']
    kcat_rga1_plus_A = results['Rga1_plus_A']
    print(f"When Rga1 is added, the kcat drops from {kcat_control}/second to {kcat_rga1}/second, indicating it is an inhibitor.")
    print(f"When excess substrate is added with Rga1, the kcat remains at {kcat_rga1_plus_A}/second.")
    if kcat_rga1_plus_A <= kcat_rga1:
        print("Conclusion: Since excess substrate does not reverse the inhibition, Rga1 is a non-competitive or irreversible inhibitor.\n")
    else:
        print("Conclusion: Rga1 is a reversible, competitive inhibitor.\n")

    print("--- Evaluating Answer Choice C ---")
    print("Choice C states: Al1 and Al2 function as allosteric modulators for Zma1. Al1 and Al2 bind to the same site on Zma1. Rga1 is an irreversible inhibitor.")

    # Part 1: Al1 and Al2 as allosteric modulators
    kcat_al2 = results['Al2']
    print(f"1. Al1 increases kcat ({kcat_control} to {kcat_al1}), so it's an activator. Al2 decreases kcat ({kcat_control} to {kcat_al2}), so it's an inhibitor. This is consistent with them being allosteric modulators. (Correct)")

    # Part 2: Al1 and Al2 binding site
    kcat_al1_al2 = results['Al1_Al2']
    print(f"2. When both Al1 and Al2 are added, the kcat is {kcat_al1_al2}/second. This rate is identical to the rate with Al2 alone ({kcat_al2}/second). This strongly suggests they compete for the same site. (Correct)")

    # Part 3: Rga1 as irreversible inhibitor
    print(f"3. As determined above, Rga1's inhibition is not overcome by substrate ({kcat_rga1}/s vs {kcat_rga1_plus_A}/s), which fits the definition of an irreversible or non-competitive inhibitor. (Correct)")

    print("\nFinal Decision: All parts of answer choice C are supported by the data.")
    print("<<<C>>>")

# Run the analysis
analyze_enzyme_data()