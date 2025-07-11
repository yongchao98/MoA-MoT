def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    (A) driving a toroidal current
    (B) elongating flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar

    The classification is as follows:
    - Only B: Axisymmetric, current-driven devices where the current (A) is the
      method to achieve the twist (B). They are planar (no C).
    - B and C: Non-axisymmetric stellarators that use 3D shaping and a
      non-planar axis (C) to create the twist (B), without a major current (no A).
    - A, B, and C: Hybrid devices using both current and non-planar geometry.
      None of the listed devices fit this category.
    """

    classification = {
        "Only property B": [
            "(1) Tokamaks",
            "(5) Reversed Field Pinches"
        ],
        "Both properties B and C": [
            "(2) LHD",
            "(3) Wendelstein 7-X",
            "(4) NCSX"
        ],
        "Properties A, B, and C": [
            "None"
        ]
    }

    print("Classification of Fusion Experiments:")
    print("-" * 35)

    for category, devices in classification.items():
        print(f"\nCategory: {category}")
        for device in devices:
            print(f"- {device}")

if __name__ == "__main__":
    classify_fusion_experiments()