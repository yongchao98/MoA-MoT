def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    Properties:
    A: Driving a toroidal current
    B: Elongating the flux surfaces and making them rotate poloidally
    C: Making the magnetic axis non-planar
    """

    # Data structure: {Experiment: [Properties]}
    experiments = {
        "Tokamaks": ['A', 'B'],
        "LHD": ['B', 'C'],
        "Wendelstein 7-X": ['B', 'C'],
        "NCSX": ['A', 'B', 'C'],
        "Reversed Field Pinches": ['A', 'B']
    }

    # Categories to find, as requested by the user
    # Note: Properties are sorted to ensure correct matching, e.g., ['C', 'B'] matches ['B', 'C']
    categories = {
        "Only property B": ['B'],
        "Properties B and C": ['B', 'C'],
        "Properties A, B, and C": ['A', 'B', 'C']
    }

    # Dictionary to hold the results
    results = {name: [] for name in categories.keys()}

    # Classify each experiment
    for exp_name, exp_props in experiments.items():
        sorted_props = sorted(exp_props)
        for cat_name, cat_props in categories.items():
            if sorted_props == sorted(cat_props):
                results[cat_name].append(exp_name)

    # Print the results
    print("Classification of Fusion Experiments:\n")
    for category, devices in results.items():
        print(f"{category}:")
        if devices:
            for device in devices:
                print(f"- {device}")
        else:
            print("- None")
        print() # Add a blank line for readability

if __name__ == "__main__":
    classify_fusion_experiments()