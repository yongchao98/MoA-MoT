def classify_fusion_experiments():
    """
    Classifies fusion experiments based on three properties for twisting the magnetic field.

    Properties:
    A: Driving a toroidal current.
    B: Elongating flux surfaces AND making them rotate poloidally around the torus.
    C: Making the magnetic axis non-planar.
    """

    # Database of experiment properties [A, B, C]
    # True if the experiment uses the property, False otherwise.
    experiments = {
        "Tokamaks": [True, False, False],
        "LHD": [False, True, True],
        "Wendelstein 7-X": [False, True, True],
        "NCSX": [True, True, True],
        "Reversed Field Pinches": [True, False, False]
    }

    # Categories as requested by the user
    # Category Name: [A, B, C] boolean pattern
    categories = {
        "Uses only property B": [False, True, False],
        "Uses both B and C": [False, True, True],
        "Uses A, B, and C": [True, True, True]
    }

    # Initialize results dictionary
    results = {name: [] for name in categories}
    unclassified = []

    # Iterate through each experiment and classify it
    for exp_name, exp_props in experiments.items():
        matched = False
        for cat_name, cat_props in categories.items():
            if exp_props == cat_props:
                results[cat_name].append(exp_name)
                matched = True
                break
        if not matched:
            unclassified.append(exp_name)

    # Print the results
    print("Classification of Fusion Experiments:")
    for cat_name, exp_list in results.items():
        print(f"{cat_name}: {exp_list}")
    
    if unclassified:
        print(f"Does not fit the specified categories: {unclassified}")


if __name__ == '__main__':
    classify_fusion_experiments()
