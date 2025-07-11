def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting magnetic fields.
    """

    # Properties based on the prompt:
    # A: driving a toroidal current
    # B: elongating flux surfaces and making them rotate poloidally
    # C: making the magnetic axis non-planar

    # Our analysis of each experiment's properties
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"A", "B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # The categories requested in the prompt
    categories = {
        "Only B": {"B"},
        "B and C": {"B", "C"},
        "A, B, and C": {"A", "B", "C"}
    }

    # Prepare lists to hold the classified experiments
    classified_results = {name: [] for name in categories}
    uncategorized = {}

    # Classify each experiment
    for exp_name, exp_props in experiments.items():
        matched = False
        for cat_name, cat_props in categories.items():
            if exp_props == cat_props:
                classified_results[cat_name].append(exp_name)
                matched = True
                break
        if not matched:
            uncategorized[exp_name] = sorted(list(exp_props))
    
    # Print the results
    print("Classification of Fusion Experiments:\n")

    for cat_name, exp_list in classified_results.items():
        print(f"Category '{cat_name}':")
        if exp_list:
            for exp in exp_list:
                print(f"  - {exp}")
        else:
            print("  - None")
        print()

    if uncategorized:
        print("The following experiments do not fit the given categories:")
        for exp_name, props in uncategorized.items():
            print(f"  - {exp_name} (Properties: {', '.join(props)})")

classify_fusion_experiments()