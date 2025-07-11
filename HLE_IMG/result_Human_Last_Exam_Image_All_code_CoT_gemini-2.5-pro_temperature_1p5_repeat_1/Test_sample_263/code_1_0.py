def analyze_emitter_stability():
    """
    Analyzes and compares the expected stability of three Ir(III) complexes
    for use in Light-Emitting Electrochemical Cells (LECs).
    """

    # Define the key stability-enhancing features for each complex.
    # We look for two main features:
    # 1. Fluorination: Increases electrochemical stability by making the complex harder to oxidize.
    # 2. Bulky Groups: Increase morphological stability by preventing aggregation.
    complex_features = {
        "Complex 1": [],
        "Complex 2": ["Bulky Groups (phenyl/tolyl)"],
        "Complex 3": ["Fluorination", "Bulky Groups (tert-butyl)"]
    }

    stability_scores = {}

    print("LEC Stability Analysis based on Emitter Molecular Structure")
    print("==========================================================")
    print("Plan: Assign a stability score to each complex. The score starts at a baseline of 1 and increases based on the presence of known stability-enhancing features.")
    print(" - Bonus for Fluorination: +2 points (enhances electrochemical stability)")
    print(" - Bonus for Bulky Groups: +2 points (enhances morphological stability)")
    print("-" * 58, "\n")

    for name, features in complex_features.items():
        score = 1  # Baseline score for a standard Ir(III) complex
        equation_parts = ["1 (baseline)"]

        if "Fluorination" in features:
            score += 2
            equation_parts.append("+ 2 (fluorination)")
        if any("Bulky Groups" in f for f in features):
            score += 2
            equation_parts.append("+ 2 (bulky groups)")

        stability_scores[name] = score
        
        print(f"--- {name} ---")
        if not features:
            print("Features: Standard benchmark structure with no special stabilizing modifications.")
        else:
            print(f"Features: {', '.join(features)}")
        
        # Outputting each number in the final equation as requested
        final_equation = " ".join(equation_parts)
        print(f"Stability Score Equation: {final_equation} = {score}")
        print("-" * (len(name) + 8), "\n")

    # Determine the most stable complex
    best_complex = max(stability_scores, key=stability_scores.get)

    print("======================\nFinal Conclusion:")
    print("Complex 3 has the highest stability score because it incorporates two distinct and effective strategies:")
    print("1. Fluorination of the cyclometalating ligand to increase resistance to oxidation.")
    print("2. Addition of bulky tert-butyl groups to the ancillary ligand to improve morphological stability.")
    print(f"\nTherefore, LECs based on {best_complex} are expected to be the most stable.")

analyze_emitter_stability()
<<<C>>>