def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting properties.

    The properties are:
    (A) driving a toroidal current
    (B) elongating flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar

    The experiments are classified into three groups:
    - Uses only property B
    - Uses both B and C
    - Uses A, B, and C
    """

    # Based on the prompt's constraints, we categorize the experiments.
    # Any device with property A (current-driven) must go in the "A, B, C" group.
    # Stellarators (no net current) are categorized based on B and C.
    categorization = {
        "Uses only property B": [
            "(None of the listed experiments fit this category)"
        ],
        "Uses both B and C": [
            "LHD",
            "Wendelstein 7-X",
            "NCSX"
        ],
        "Uses A, B, and C": [
            "Tokamaks",
            "Reversed Field Pinches"
        ]
    }

    print("Categorization of Fusion Experiments:")
    print("-" * 40)

    for category, experiments in categorization.items():
        print(f"{category}:")
        for experiment in experiments:
            print(f"  - {experiment}")
        print("") # Add a newline for better readability

if __name__ == '__main__':
    classify_fusion_experiments()