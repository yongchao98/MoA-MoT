def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their methods for twisting magnetic fields.

    The properties are:
    A) driving a toroidal current
    B) elongating the flux surfaces and making them rotate poloidally
    C) making the magnetic axis non-planar
    """

    # Dictionary mapping each experiment to its classification.
    # The numbering corresponds to the user's request.
    classifications = {
        "1) Tokamaks": "Only B",
        "2) LHD": "Both B and C",
        "3) Wendelstein 7-X": "Both B and C",
        "4) NCSX": "Both B and C",
        "5) Reversed Field Pinches": "A, B, and C"
    }

    print("Fusion Experiment Classifications:")
    print("-" * 35)

    # Print the classification for each numbered experiment.
    for experiment, category in classifications.items():
        print(f"{experiment}: Uses property/properties {category}")

    print("-" * 35)
    print("\nExplanation:")
    print(" - Tokamak: While driven by a current (A), the resulting twisted field structure is property (B). In this simplified scheme and with a planar axis, it's categorized as 'Only B'.")
    print(" - LHD, W7-X, NCSX: These are all advanced stellarators that use complex external coils to create shaped fields (B) on a non-planar axis (C).")
    print(" - Reversed Field Pinches: Use a strong current (A) creating a twisted field (B), and rely on an inherent 3D dynamo effect, which acts as a 'non-planar' feature (C).")

if __name__ == "__main__":
    classify_fusion_experiments()
    # The final answer can be represented by concatenating the category labels for each device 1 through 5.
    final_answer_string = "Only B; Both B and C; Both B and C; Both B and C; A, B, and C"
    # To conform to the specified output format, we must provide a single answer string.
    # The requested format "<<<answer content>>>" will contain the classification for all five devices.
    print(f"\n<<<Only B, Both B and C, Both B and C, Both B and C, A, B, and C>>>")
