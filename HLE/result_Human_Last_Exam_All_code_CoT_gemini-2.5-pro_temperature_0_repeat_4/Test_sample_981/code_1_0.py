def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field properties
    and prints the results.
    """
    # Define the properties for each fusion experiment:
    # A: Driving a toroidal current
    # B: Elongating flux surfaces and poloidal rotation (twist)
    # C: Non-planar magnetic axis
    experiments = {
        "Tokamaks": {"A": True, "B": True, "C": False},
        "LHD": {"A": False, "B": True, "C": True},
        "Wendelstein 7-X": {"A": False, "B": True, "C": True},
        "NCSX": {"A": False, "B": True, "C": True},
        "Reversed Field Pinches": {"A": True, "B": True, "C": False}
    }

    # The categories from the prompt are ("Only B", "B and C", "A, B, and C").
    # This set is inconsistent with the physics of the devices listed.
    # We assume a typo and use a corrected set of categories for a full classification.
    # Specifically, we assume "A, B, and C" was meant to be "A and B".
    
    classification = {
        "Uses properties B and C": [],
        "Uses properties A and B": [],
        "Uses only property B": [],
        "Uses properties A, B, and C": []
    }

    # Classify each experiment based on its properties
    for name, props in experiments.items():
        if props["A"] and props["B"] and props["C"]:
            classification["Uses properties A, B, and C"].append(name)
        elif props["B"] and props["C"]:
            classification["Uses properties B and C"].append(name)
        elif props["A"] and props["B"]:
            classification["Uses properties A and B"].append(name)
        elif props["B"]:
            classification["Uses only property B"].append(name)

    # Print the results in a clear format
    print("Classification of Fusion Experiments:")
    print("="*40)

    # Category: B and C
    print("\nExperiments that use both B and C:")
    print("(Elongated/rotating surfaces and a non-planar axis)")
    if classification["Uses properties B and C"]:
        for device in classification["Uses properties B and C"]:
            print(f"- {device}")
    else:
        print("- None")

    # Category: A and B (The logical category for Tokamaks/RFPs)
    print("\nExperiments that use both A and B:")
    print("(Toroidal current and elongated/rotating surfaces)")
    if classification["Uses properties A and B"]:
        for device in classification["Uses properties A and B"]:
            print(f"- {device}")
    else:
        print("- None")
        
    # Category: Only B (As requested by prompt)
    print("\nExperiments that use only property B:")
    if classification["Uses only property B"]:
        for device in classification["Uses only property B"]:
            print(f"- {device}")
    else:
        print("- None")

    # Category: A, B, and C (As requested by prompt)
    print("\nExperiments that use A, B, and C:")
    if classification["Uses properties A, B, and C"]:
        for device in classification["Uses properties A, B, and C"]:
            print(f"- {device}")
    else:
        print("- None")
        
    print("\n" + "="*40)
    print("Note: The prompt's categories ('Only B', 'B and C', 'A, B, and C') do not fully cover the listed devices.")
    print("Tokamaks and RFPs use properties A and B. The classification above shows where each device fits based on its physical principles.")

classify_fusion_experiments()