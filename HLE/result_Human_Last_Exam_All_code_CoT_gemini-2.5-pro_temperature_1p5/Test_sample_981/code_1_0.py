def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting properties.
    A: Driving a toroidal current
    B: Elongating flux surfaces with poloidal rotation
    C: Non-planar magnetic axis
    """
    # Define the properties for each experiment using sets for easy comparison.
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"A", "B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # Define the target property sets for classification.
    cat_only_b = {"B"}
    cat_b_and_c = {"B", "C"}
    cat_a_b_and_c = {"A", "B", "C"}

    # Initialize lists to hold the names for each category.
    only_b_list = []
    b_and_c_list = []
    a_b_and_c_list = []

    # Classify each experiment by comparing its properties to the target sets.
    for name, properties in experiments.items():
        if properties == cat_only_b:
            only_b_list.append(name)
        elif properties == cat_b_and_c:
            b_and_c_list.append(name)
        elif properties == cat_a_b_and_c:
            a_b_and_c_list.append(name)

    # Print the results in a clear format.
    print("Classification of Fusion Experiments:\n")

    print("1. Experiments that use ONLY property B:")
    if only_b_list:
        for exp in only_b_list:
            print(f"- {exp}")
    else:
        print("- None")

    print("\n2. Experiments that use BOTH properties B and C:")
    if b_and_c_list:
        for exp in b_and_c_list:
            print(f"- {exp}")
    else:
        print("- None")

    print("\n3. Experiments that use properties A, B, AND C:")
    if a_b_and_c_list:
        for exp in a_b_and_c_list:
            print(f"- {exp}")
    else:
        print("- None")

if __name__ == "__main__":
    classify_fusion_experiments()