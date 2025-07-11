def classify_fusion_experiments():
    """
    Classifies fusion experiments based on how they twist their magnetic fields,
    following a specific interpretation to fit them into the given categories.
    """

    print("Classifying fusion experiments based on the properties used to twist the magnetic field.")
    print("The properties are (A) toroidal current, (B) external coil shaping, and (C) non-planar axis.")
    print("-" * 20)

    # The list of experiments as given in the problem
    experiments = {
        1: "Tokamaks",
        2: "LHD",
        3: "Wendelstein 7-X",
        4: "NCSX",
        5: "Reversed Field Pinches"
    }

    # Apply the classification logic developed in the plan.
    only_B_names = ["LHD"]
    B_and_C_names = ["Wendelstein 7-X", "NCSX"]
    A_B_and_C_names = ["Tokamaks", "Reversed Field Pinches"]

    # Create lists of numbered experiment strings for the output.
    only_B = [f"{num}) {name}" for num, name in experiments.items() if name in only_B_names]
    B_and_C = [f"{num}) {name}" for num, name in experiments.items() if name in B_and_C_names]
    A_B_and_C = [f"{num}) {name}" for num, name in experiments.items() if name in A_B_and_C_names]

    print("Final Classification:\n")

    print("Use only property B (Stellarators whose primary feature is helical windings):")
    print(f"  {', '.join(only_B)}\n")

    print("Use both B and C (Stellarators primarily designed with an optimized non-planar axis):")
    print(f"  {', '.join(B_and_C)}\n")
    
    print("Use A, B, and C (Current-driven devices):")
    print(f"  {', '.join(A_B_and_C)}")

# Execute the function to print the classification.
classify_fusion_experiments()