def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their methods for twisting magnetic fields.

    Properties:
    A: Driving a toroidal current.
    B: Elongating the flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.
    """
    # Step 1 & 2: Define experiment properties in a dictionary
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"A", "B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # Define the target categories as sets of properties
    categories = {
        "only property B": {"B"},
        "both B and C": {"B", "C"},
        "A, B, and C": {"A", "B", "C"}
    }

    # Step 3 & 4: Categorize the experiments
    results = {key: [] for key in categories.keys()}

    for exp_name, exp_props in experiments.items():
        for cat_name, cat_props in categories.items():
            if exp_props == cat_props:
                results[cat_name].append(exp_name)

    # Step 5: Print the formatted results
    print("Classification of fusion experiments based on magnetic field properties:\n")
    for category, exp_list in results.items():
        # Format the list of experiments for printing
        if not exp_list:
            experiment_str = "None"
        else:
            experiment_str = ", ".join(exp_list)
        
        print(f"Experiments that use {category}: {experiment_str}")

# Execute the classification
classify_fusion_experiments()

# Final Answer Summary
# Use only property B: None
# Use both B and C: LHD, Wendelstein 7-X
# Use A, B, and C: NCSX
# The other experiments (Tokamaks, Reversed Field Pinches) use properties A and B.
# <<<Use only B: None; Use both B and C: LHD, Wendelstein 7-X; Use A, B, and C: NCSX>>>