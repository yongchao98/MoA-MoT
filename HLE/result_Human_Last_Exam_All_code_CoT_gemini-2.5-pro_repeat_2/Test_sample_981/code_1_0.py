def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    (A) driving a toroidal current
    (B) elongating the flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar

    The allowed classifications are:
    - Only B
    - B and C
    - A, B, and C
    """

    # This dictionary maps each experiment's number and name to its classification.
    experiments = {
        1: {"name": "Tokamaks", "category_key": "Only B"},
        2: {"name": "LHD", "category_key": "B and C"},
        3: {"name": "Wendelstein 7-X", "category_key": "B and C"},
        4: {"name": "NCSX", "category_key": "A, B, and C"},
        5: {"name": "Reversed Field Pinches", "category_key": "Only B"},
    }

    # Group experiments by their classification category.
    classifications = {
        "Only B": [],
        "B and C": [],
        "A, B, and C": [],
    }

    for number, data in experiments.items():
        classifications[data["category_key"]].append(f"({number}) {data['name']}")

    # Print the results in a structured format.
    print("Classification of Fusion Experiments:")
    print("-" * 35)

    print("\nExperiments using only property B:")
    for exp in classifications["Only B"]:
        print(f"- {exp}")

    print("\nExperiments using both properties B and C:")
    for exp in classifications["B and C"]:
        print(f"- {exp}")

    print("\nExperiments using properties A, B, and C:")
    for exp in classifications["A, B, and C"]:
        print(f"- {exp}")

if __name__ == "__main__":
    classify_fusion_experiments()