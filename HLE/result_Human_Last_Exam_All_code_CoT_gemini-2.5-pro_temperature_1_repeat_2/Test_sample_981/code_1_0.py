def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    Properties:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally (shaping/twist).
    C: Making the magnetic axis non-planar.
    """

    # Step 1: Define the properties for each experiment
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # Step 2: Initialize lists for the categories
    # Note: As analyzed, no device fits "A, B, and C". We assume the category
    # for Tokamaks/RFPs is "A and B". No listed device fits "only B".
    only_B = []
    B_and_C = []
    A_and_B = [] # Corrected category based on analysis
    A_B_and_C = []

    # Step 3: Iterate and classify experiments
    for exp, properties in experiments.items():
        if properties == {"B"}:
            only_B.append(exp)
        elif properties == {"B", "C"}:
            B_and_C.append(exp)
        elif properties == {"A", "B"}:
            A_and_B.append(exp)
        elif properties == {"A", "B", "C"}:
            A_B_and_C.append(exp)

    # Step 4: Print the results in a structured format
    print("Classification of Fusion Experiments:\n")

    print("Uses both B and C (Shaping/Rotation from coils and Non-planar Axis):")
    if B_and_C:
        for exp in B_and_C:
            print(f"- {exp}")
    else:
        print("- None in this list")
    
    # This category is created to house the remaining experiments, as it's the
    # most logical fit not explicitly mentioned in the prompt's categories.
    print("\nUses both A and B (Toroidal Current and Shaping/Rotation):")
    if A_and_B:
        for exp in A_and_B:
            print(f"- {exp}")
    else:
        print("- None in this list")

    print("\nUses only B (Shaping/Rotation from coils, Planar Axis):")
    if only_B:
        for exp in only_B:
            print(f"- {exp}")
    else:
        print("- None in this list")
    
    print("\nUses A, B, and C (Current, Shaping/Rotation, and Non-planar Axis):")
    if A_B_and_C:
        for exp in A_B_and_C:
            print(f"- {exp}")
    else:
        print("- None in this list")


if __name__ == "__main__":
    classify_fusion_experiments()