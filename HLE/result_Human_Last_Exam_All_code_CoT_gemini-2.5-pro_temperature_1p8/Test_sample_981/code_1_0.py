def classify_fusion_experiments():
    """
    Classifies fusion experiments based on the method used to twist the magnetic field.

    The properties are:
    (A) Driving a toroidal current.
    (B) Elongating flux surfaces and making them rotate poloidally.
    (C) Making the magnetic axis non-planar.
    """

    # Classification based on the analysis
    # B and C: Stellarators, which use 3D shaping and a non-planar axis.
    # A, B, and C: Current-driven devices, which also use shaping and can have 3D field components.
    
    only_B = []
    uses_B_and_C = ["LHD", "Wendelstein 7-X", "NCSX"]
    uses_A_B_and_C = ["Tokamaks", "Reversed Field Pinches"]

    print("Classification of Fusion Experiments by Magnetic Field Twisting Method:")
    print("=" * 65)

    print("\nCategory: Uses only property B (shaping/rotation on a planar axis, no current)")
    if not only_B:
        print("  - None of the listed experiments fit this category.")
    else:
        for experiment in only_B:
            print(f"  - {experiment}")

    print("\nCategory: Uses both properties B and C (shaping/rotation and a non-planar axis)")
    if not uses_B_and_C:
        print("  - None of the listed experiments fit this category.")
    else:
        for experiment in uses_B_and_C:
            print(f"  - {experiment}")

    print("\nCategory: Uses properties A, B, and C (current, shaping/rotation, and non-planar features)")
    if not uses_A_B_and_C:
        print("  - None of the listed experiments fit this category.")
    else:
        for experiment in uses_A_B_and_C:
            print(f"  - {experiment}")

    print("\n" + "=" * 65)

# Execute the classification and print the results
classify_fusion_experiments()

# The following is a placeholder for the final answer format as requested.
# The textual output from the code above is the primary answer.
final_answer = {
    "Uses only property B": [],
    "Uses both B and C": ["LHD", "Wendelstein 7-X", "NCSX"],
    "Uses A, B, and C": ["Tokamaks", "Reversed Field Pinches"]
}
# <<<final_answer>>> is a placeholder and should not be printed literally.
# The code's printout serves as the answer.