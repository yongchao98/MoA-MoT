def classify_fusion_experiments():
    """
    Classifies fusion experiments based on how they twist their magnetic fields
    according to the provided properties and categories.
    """

    # Properties:
    # A: driving a toroidal current
    # B: elongating flux surfaces and making them rotate poloidally
    # C: making the magnetic axis non-planar
    
    # Categories for classification
    only_b = []
    b_and_c = []
    a_b_and_c = []

    # Full list of experiments
    experiments = {
        "Tokamaks": "Only B",
        "LHD": "B and C",
        "Wendelstein 7-X": "B and C",
        "NCSX": "B and C",
        "Reversed Field Pinches": "A, B, and C"
    }

    # Sort experiments into their categories
    for exp, category in experiments.items():
        if category == "Only B":
            only_b.append(exp)
        elif category == "B and C":
            b_and_c.append(exp)
        elif category == "A, B, and C":
            a_b_and_c.append(exp)

    # Print the final classification
    print("Classification of Fusion Experiments:")
    print("-" * 35)
    
    print("Use Only Property B:")
    for exp in only_b:
        print(f"- {exp}")
    
    print("\nUse Both Properties B and C:")
    for exp in b_and_c:
        print(f"- {exp}")

    print("\nUse Properties A, B, and C:")
    for exp in a_b_and_c:
        print(f"- {exp}")

if __name__ == "__main__":
    classify_fusion_experiments()