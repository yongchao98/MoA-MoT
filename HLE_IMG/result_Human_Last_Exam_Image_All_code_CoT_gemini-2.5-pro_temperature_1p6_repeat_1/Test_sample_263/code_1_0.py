def analyze_lec_stability():
    """
    Analyzes the expected stability of LECs based on three different Ir(III) complexes.
    Stability is scored based on key molecular features that prevent degradation.
    """

    # --- Step 1: Define scoring criteria for stability ---
    # Score points are assigned for features known to enhance stability in emitters.
    # - Bulky groups prevent aggregation (morphological stability).
    # - Electron-Withdrawing Groups (EWGs) prevent oxidative degradation (electrochemical stability).
    scoring_rules = {
        'bulky_groups': 1,
        'electron_withdrawing_groups': 2  # Often has a more significant impact on lifetime
    }

    # --- Step 2: Define the properties of each complex based on the image ---
    complex_data = {
        'Complex 1': {'has_bulky_groups': False, 'has_ewgs': False},
        'Complex 2': {'has_bulky_groups': True, 'has_ewgs': False},
        'Complex 3': {'has_bulky_groups': True, 'has_ewgs': True}
    }

    print("--- Analysis of Emitter Stability for LECs ---")
    print("Scoring based on key stabilizing features:")
    print(f" - Bulky Groups: +{scoring_rules['bulky_groups']} pt")
    print(f" - Electron-Withdrawing Groups (EWGs): +{scoring_rules['electron_withdrawing_groups']} pts\n")

    # --- Step 3: Calculate and display the stability score for each complex ---
    results = {}
    for name, features in complex_data.items():
        score = 0
        calculation_str_parts = []
        
        if features['has_bulky_groups']:
            score += scoring_rules['bulky_groups']
            calculation_str_parts.append(str(scoring_rules['bulky_groups']))
        
        if features['has_ewgs']:
            score += scoring_rules['electron_withdrawing_groups']
            calculation_str_parts.append(str(scoring_rules['electron_withdrawing_groups']))
            
        # To handle the case with a score of 0
        if not calculation_str_parts:
            calculation_str = "0"
        else:
            calculation_str = " + ".join(calculation_str_parts)

        results[name] = score
        
        print(f"Analysis for {name}:")
        print(f"  - Features: Bulky Groups={features['has_bulky_groups']}, EWGs={features['has_ewgs']}")
        print(f"  - Stability Score Equation: {calculation_str} = {score}\n")

    # --- Step 4: Determine the most stable complex ---
    most_stable_complex = max(results, key=results.get)
    
    print("--- Conclusion ---")
    print(f"The most stable complex is {most_stable_complex}, with the highest stability score.")
    print("It benefits from both steric hindrance (bulky groups) and enhanced electrochemical stability (EWGs).")
    print("Therefore, LECs based on Complex 3 are expected to be the most stable.")

# Run the analysis
analyze_lec_stability()