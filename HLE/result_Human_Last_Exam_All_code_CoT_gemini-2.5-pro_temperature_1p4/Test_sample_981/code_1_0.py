def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their methods for twisting the magnetic field.

    Properties:
    A: Driving a toroidal current.
    B: Elongating flux surfaces and making them rotate poloidally.
    C: Making the magnetic axis non-planar.
    """

    # Data structure: {Category: [List of experiments]}
    classification = {
        "Only B": [],
        "B and C": ["LHD", "Wendelstein 7-X"],
        "A, B, and C": ["NCSX"]
    }

    # Although not requested by the prompt, this is where other devices would go
    # for completeness in the logic.
    # classification["A and B"] = ["Tokamaks", "Reversed Field Pinches"]

    print("Classification of Fusion Experiments:\n")

    category_b_only = "Uses only property B:"
    print(category_b_only)
    if not classification["Only B"]:
        print(" - None")
    else:
        for experiment in classification["Only B"]:
            print(f" - {experiment}")
    print("-" * 30)

    category_b_and_c = "Uses both B and C:"
    print(category_b_and_c)
    if not classification["B and C"]:
        print(" - None")
    else:
        for experiment in classification["B and C"]:
            print(f" - {experiment}")
    print("-" * 30)

    category_a_b_and_c = "Uses A, B, and C:"
    print(category_a_b_and_c)
    if not classification["A, B, and C"]:
        print(" - None")
    else:
        for experiment in classification["A, B, and C"]:
            print(f" - {experiment}")
    print("-" * 30)


if __name__ == "__main__":
    classify_fusion_experiments()
    # The final answer is the conceptual classification derived from physics principles.
    # The code serves to present this classification clearly.
    # The user is expected to understand the categories from the printed output.
    # To conform to the output format, I will serialize the classification.
    answer_string = "Only B: None; B and C: LHD, Wendelstein 7-X; A, B, and C: NCSX"
    # print(f"<<<{answer_string}>>>")