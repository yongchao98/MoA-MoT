def analyze_complex_stability():
    """
    Analyzes and compares the expected stability of three Ir(III) complexes
    for use in Light-emitting Electrochemical Cells (LECs).
    """
    # Define complexes with flags for key stability-enhancing features.
    # Features: 'fluorinated' (strengthens Ir-C bond), 'bulky_groups' (provides steric protection)
    complexes = {
        "Complex 1": {"fluorinated": False, "bulky_groups": False},
        "Complex 2": {"fluorinated": False, "bulky_groups": False}, # dppz ligand is large but not sterically bulky in a protective way
        "Complex 3": {"fluorinated": True, "bulky_groups": True}
    }

    print("Evaluating LEC emitter stability based on molecular structure.")
    print("A scoring system is used to quantify stability enhancements:")
    print("- Fluorination of cyclometalating ligand: +2 points")
    print("- Bulky groups on ancillary ligand: +2 points")
    print("-" * 50)

    results = {}
    for name, features in complexes.items():
        score = 0
        
        # Score contribution from fluorination
        fluorination_score = 2 if features["fluorinated"] else 0
        score += fluorination_score
        
        # Score contribution from bulky groups
        bulky_group_score = 2 if features["bulky_groups"] else 0
        score += bulky_group_score
        
        results[name] = {
            "score": score,
            "fluorination_score": fluorination_score,
            "bulky_group_score": bulky_group_score
        }

    # Print the analysis for each complex
    for name, data in results.items():
        print(f"Analysis for {name}:")
        print(f"Stability Score = {data['fluorination_score']} (from fluorination) + {data['bulky_group_score']} (from bulky groups) = {data['score']}")
        print("")

    # Determine the most stable complex
    most_stable_complex = max(results, key=lambda k: results[k]['score'])
    
    print("-" * 50)
    print("Conclusion:")
    print(f"{most_stable_complex} has the highest stability score.")
    print("The combination of fluorination and bulky groups provides superior electrochemical and morphological stability.")
    print("Therefore, LECs based on Complex 3 are expected to be the most stable.")

# Run the analysis
analyze_complex_stability()