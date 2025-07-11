def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field properties.

    Properties:
    A: Driving or possessing a significant toroidal current.
    B: Elongating flux surfaces and creating poloidal rotation (rotational transform).
    C: Having a non-planar magnetic axis.
    """

    # Data structure mapping each experiment to its set of properties.
    # Note on Property A: This includes significant intrinsic bootstrap current,
    # which is why LHD and NCSX are included, while W7-X (optimized to
    # minimize this current) is not.
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"A", "B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"A", "B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # The categories requested by the user
    target_categories = {
        "Only B": {"B"},
        "Both B and C": {"B", "C"},
        "A, B, and C": {"A", "B", "C"}
    }

    # A dictionary to hold the final classification
    classification = {key: [] for key in target_categories}

    # Classify each experiment
    for name, properties in experiments.items():
        for category_name, category_props in target_categories.items():
            if properties == category_props:
                classification[category_name].append(name)
                break

    # Print the results in a formatted way
    print("Classification of Fusion Experiments by Magnetic Field Properties:")
    print("-" * 60)
    print("(A) Toroidal Current, (B) Shaping/Rotation, (C) Non-planar Axis\n")

    category_descriptions = {
        "Only B": "Use only property B:",
        "Both B and C": "Use both properties B and C:",
        "A, B, and C": "Use properties A, B, and C:"
    }
    
    for category_name, description in category_descriptions.items():
        print(description)
        results = classification[category_name]
        if results:
            for experiment in sorted(results):
                print(f"- {experiment}")
        else:
            print("- None of the listed experiments fall into this category.")
        print()

# Execute the classification and print the output
classify_fusion_experiments()