def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    Properties:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.
    """
    experiments = [
        {'name': 'Tokamaks', 'A': True, 'B': True, 'C': False},
        {'name': 'LHD', 'A': False, 'B': True, 'C': True},
        {'name': 'Wendelstein 7-X', 'A': False, 'B': True, 'C': True},
        {'name': 'NCSX', 'A': True, 'B': True, 'C': True},
        {'name': 'Reversed Field Pinches', 'A': True, 'B': True, 'C': False}
    ]

    # Categories based on the prompt
    categories = {
        "Only B": [],
        "B and C": [],
        "A, B, and C": [],
        "Uncategorized": []
    }

    for exp in experiments:
        props = (exp['A'], exp['B'], exp['C'])
        if not exp['A'] and exp['B'] and not exp['C']:
            categories["Only B"].append(exp['name'])
        elif not exp['A'] and exp['B'] and exp['C']:
            categories["B and C"].append(exp['name'])
        elif exp['A'] and exp['B'] and exp['C']:
            categories["A, B, and C"].append(exp['name'])
        else:
            # Store properties for uncategorized for clarity
            uncategorized_props = []
            if exp['A']: uncategorized_props.append('A')
            if exp['B']: uncategorized_props.append('B')
            if exp['C']: uncategorized_props.append('C')
            categories["Uncategorized"].append(f"{exp['name']} (Properties: {', '.join(uncategorized_props)})")

    print("Classification of Fusion Experiments:")
    for category, members in categories.items():
        # Adjusting the label for the "Uncategorized" group for better readability
        if category == "Uncategorized":
            print("\nExperiments that do not fit the given categories:")
        else:
            print(f"\nCategory '{category}':")

        if members:
            for member in members:
                print(f"- {member}")
        else:
            print("- None")

if __name__ == "__main__":
    classify_fusion_experiments()