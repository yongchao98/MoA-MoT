def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.
    A: driving a toroidal current
    B: elongating the flux surfaces and making them rotate poloidally
    C: making the magnetic axis non-planar
    """

    # Define the properties for each experiment based on physics principles.
    # 'A' = Toroidal Current, 'B' = Shaping/Rotation, 'C' = Non-planar Axis
    experiments = {
        'Tokamaks': {'A', 'B'},
        'LHD': {'B', 'C'},
        'Wendelstein 7-X': {'B', 'C'},
        'NCSX': {'A', 'B', 'C'},
        'Reversed Field Pinches': {'A', 'B'}
    }

    # Define the categories provided in the problem description.
    categories = {
        'Uses only property B': {'B'},
        'Uses both B and C': {'B', 'C'},
        'Uses A, B, and C': {'A', 'B', 'C'}
    }

    # Prepare a dictionary to hold the classification results.
    # Initialize with the defined categories and an extra one for non-fitting items.
    results = {name: [] for name in categories.keys()}
    uncategorized_group_name = "Uncategorized (Use properties A and B, which is not a listed category)"
    results[uncategorized_group_name] = []


    # Iterate through each experiment and assign it to a category.
    for device, properties in experiments.items():
        found_category = False
        for cat_name, cat_properties in categories.items():
            if properties == cat_properties:
                results[cat_name].append(device)
                found_category = True
                break
        if not found_category:
            results[uncategorized_group_name].append(device)

    # Print the final classification.
    print("Classification of Fusion Experiments:\n")
    for category, devices in results.items():
        # Only print categories that have devices.
        if devices:
            print(f"{category}:")
            for device in devices:
                print(f"- {device}")
            print("-" * 20)

# Execute the classification.
classify_fusion_experiments()

# The final answer summary based on the script's output.
final_answer = "Uses only property B: None. Uses both B and C: LHD, Wendelstein 7-X. Uses A, B, and C: NCSX. Uncategorized (A and B): Tokamaks, Reversed Field Pinches."
# print(f"\n<<<ANSER>>>\n{final_answer}\n<<<ANSER>>>")