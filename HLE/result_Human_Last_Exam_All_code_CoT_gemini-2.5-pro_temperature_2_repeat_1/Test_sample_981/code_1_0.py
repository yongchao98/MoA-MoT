def classify_fusion_experiments():
    """
    Classifies fusion experiments based on how they twist their magnetic fields.

    Properties:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.
    """

    # Categorization based on the analysis
    only_B = []
    B_and_C = ["LHD", "Wendelstein 7-X", "NCSX"]
    A_B_and_C = ["Tokamaks", "Reversed Field Pinches"]

    print("Classification of Fusion Experiments:\n")

    # Print the category for 'only property B'
    print("Experiments that use only property B (elongating/rotating flux surfaces):")
    if not only_B:
        print(" - None of the listed experiments fall into this category.")
    else:
        for exp in only_B:
            print(f" - {exp}")
    print("-" * 40)

    # Print the category for 'B and C'
    print("Experiments that use both properties B and C (non-planar axis and shaping):")
    for exp in B_and_C:
        print(f" - {exp}")
    print("-" * 40)

    # Print the category for 'A, B, and C'
    print("Experiments that use properties A, B, and C (toroidal current, shaping, and non-planar elements):")
    for exp in A_B_and_C:
        print(f" - {exp}")
    print("-" * 40)

if __name__ == '__main__':
    classify_fusion_experiments()