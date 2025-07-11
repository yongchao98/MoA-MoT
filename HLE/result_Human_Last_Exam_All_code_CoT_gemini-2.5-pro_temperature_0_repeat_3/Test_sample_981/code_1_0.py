def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting properties.

    Properties:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.
    """

    # Data mapping experiments to their properties
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"A", "B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # Initialize lists for each category
    only_b = []
    b_and_c = []
    a_b_and_c = []

    # Categorize each experiment
    for name, properties in experiments.items():
        if properties == {"B"}:
            only_b.append(name)
        elif properties == {"B", "C"}:
            b_and_c.append(name)
        elif properties == {"A", "B", "C"}:
            a_b_and_c.append(name)

    # Print the results
    print("Classification of Fusion Experiments:")
    print("-" * 35)

    # Format and print the list for "Only B"
    print("Use only property B:")
    if only_b:
        print(f"  - {', '.join(only_b)}")
    else:
        print("  - None")

    # Format and print the list for "B and C"
    print("\nUse both B and C:")
    if b_and_c:
        print(f"  - {', '.join(b_and_c)}")
    else:
        print("  - None")

    # Format and print the list for "A, B, and C"
    print("\nUse A, B, and C:")
    if a_b_and_c:
        print(f"  - {', '.join(a_b_and_c)}")
    else:
        print("  - None")

if __name__ == "__main__":
    classify_fusion_experiments()
    # The final answer is derived from the output of the script.
    final_answer = "Use only property B: None\nUse both B and C: LHD, Wendelstein 7-X\nUse A, B, and C: NCSX"
    # print(f"\n<<<{final_answer}>>>") # This line is for internal processing and would not be in the final user-facing code.