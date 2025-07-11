def solve_clinical_case():
    """
    Analyzes a clinical scenario to determine the most indicative lab parameter
    for a rapid decline in renal function in a patient with suspected SLE.
    """

    # Define potential lab tests and their clinical significance in this context.
    # The 'score' is a heuristic to quantify relevance to the specific question asked.
    # We want a marker for an ACTIVE, CAUSAL process.
    lab_tests = [
        {
            "name": "Serum Creatinine/BUN",
            "significance": "Measures kidney filtration function. High levels are a RESULT of kidney damage, not the direct cause.",
            "score": 1  # Indicates a problem, but not the cause.
        },
        {
            "name": "Anti-dsDNA Antibody Titer",
            "significance": "Measures autoantibodies that DRIVE the disease process. A rising level indicates high potential for a flare.",
            "score": 3  # A primary driver, very important.
        },
        {
            "name": "Complement Levels (C3 & C4)",
            "significance": "These proteins are CONSUMED during active immune-complex deposition in the kidneys. Low levels are a direct marker of an ONGOING immunological attack causing tissue damage.",
            "score": 5  # Best reflects the *active, real-time* pathological process causing the damage.
        },
        {
            "name": "Urinalysis with Microscopy",
            "significance": "Shows evidence of glomerulonephritis (e.g., red blood cell casts). It CONFIRMS the type of damage but doesn't point to the systemic cause.",
            "score": 2  # Confirms local damage, but not the systemic cause.
        }
    ]

    # Find the test with the highest score
    best_indicator = max(lab_tests, key=lambda x: x['score'])

    print("Analysis of Clinical Scenario:")
    print("1. The patient's history suggests an underlying autoimmune disease, likely Systemic Lupus Erythematosus (SLE).")
    print("2. The rapid deterioration of kidney function after an inflammatory rebound is characteristic of a severe lupus nephritis flare.")
    print("3. The key question is to identify the lab marker that best indicates the CAUSE of this rapid decline.")
    print("\nEvaluation of Lab Parameters:")
    
    for test in sorted(lab_tests, key=lambda x: x['score'], reverse=True):
        print(f"- {test['name']}: {test['significance']}")

    print("\n-------------------------------------------------------------")
    print("Conclusion:")
    print(f"The lab parameter that could have best indicated the cause of the rapid renal function decline is: {best_indicator['name']}")
    print(f"Reasoning: {best_indicator['significance']}")
    print("-------------------------------------------------------------")

solve_clinical_case()