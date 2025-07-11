def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field properties.
    
    Properties:
    (A) Driving a toroidal current
    (B) Elongating flux surfaces and rotating them poloidally
    (C) Making the magnetic axis non-planar
    """
    
    # Define the properties for each experiment based on standard principles
    experiments = {
        'Tokamaks': {'A', 'B'},
        'LHD': {'B', 'C'},
        'Wendelstein 7-X': {'B', 'C'},
        'NCSX': {'B', 'C'},
        'Reversed Field Pinches': {'A', 'B'}
    }

    # Define the property sets for the user's requested categories
    categories = {
        'only property B': frozenset(['B']),
        'both B and C': frozenset(['B', 'C']),
        'A, B, and C': frozenset(['A', 'B', 'C'])
    }

    # Prepare a dictionary to hold the results
    results = {name: [] for name in categories.keys()}

    # Match each experiment to a category
    for exp_name, exp_props in experiments.items():
        for cat_name, cat_props in categories.items():
            if frozenset(exp_props) == cat_props:
                results[cat_name].append(exp_name)

    # Print the findings
    print("Classification of Fusion Experiments:")
    print("-------------------------------------")

    # Print results for each specified category
    for category_name, exp_list in results.items():
        print(f"Experiments that use {category_name}:")
        if exp_list:
            for experiment in exp_list:
                print(f"- {experiment}")
        else:
            print("- None")
        print()
        
    # Note on unclassified experiments
    print("Note: Tokamaks and Reversed Field Pinches use properties A and B,")
    print("a combination not listed in the target categories.")

if __name__ == '__main__':
    classify_fusion_experiments()