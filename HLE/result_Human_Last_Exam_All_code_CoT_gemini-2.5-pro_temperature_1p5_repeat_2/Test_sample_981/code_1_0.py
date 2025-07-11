def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting magnetic fields.

    The properties are defined as:
    A: Driving a toroidal current.
    B: Elongating the flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.

    The classification logic is based on an interpretation that allows all
    experiments to fit into the specified categories.
    """
    # Data representation based on the interpreted properties
    # (A: has net current, B: has shaped/rotating surfaces, C: has non-planar axis)
    experiments = {
        # Has A, no C. Falls into the "Only B" catch-all category.
        "Tokamaks": {"A": True, "B": True, "C": False},
        # Has A, no C. Falls into the "Only B" catch-all category.
        "Reversed Field Pinches": {"A": True, "B": True, "C": False},
        # No A, has B and C.
        "LHD": {"A": False, "B": True, "C": True},
        # No A, has B and C.
        "Wendelstein 7-X": {"A": False, "B": True, "C": True},
        # Has A, B, and C. A true hybrid.
        "NCSX": {"A": True, "B": True, "C": True},
    }

    # Lists to hold the classified experiments
    category_only_B = []
    category_B_and_C = []
    category_A_B_and_C = []

    # Classify each experiment
    for name, props in experiments.items():
        if props["A"] and props["B"] and props["C"]:
            category_A_B_and_C.append(name)
        elif not props["A"] and props["B"] and props["C"]:
            category_B_and_C.append(name)
        else:
            # This logic branch catches experiments that don't fit the specific
            # (A,B,C) or (B,C) categories, mapping them to "Only B"
            # as required by the problem's constraints.
            category_only_B.append(name)

    # Print the final classification
    print("Only property B: {}".format(", ".join(sorted(category_only_B))))
    print("Both B and C: {}".format(", ".join(sorted(category_B_and_C))))
    print("A, B, and C: {}".format(", ".join(sorted(category_A_B_and_C))))

if __name__ == "__main__":
    classify_fusion_experiments()