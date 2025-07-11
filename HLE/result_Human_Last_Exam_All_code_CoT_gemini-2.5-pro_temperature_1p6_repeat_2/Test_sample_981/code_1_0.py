def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method for twisting magnetic fields.

    The properties are:
    (A) driving a toroidal current
    (B) elongating flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar
    """
    
    # Based on the analysis, the experiments are grouped into the given categories.
    classification = {
        "Uses only property B": [
            "LHD"
        ],
        "Uses both B and C": [
            "Wendelstein 7-X",
            "NCSX"
        ],
        "Uses A, B, and C": [
            "Tokamaks",
            "Reversed Field Pinches"
        ]
    }

    print("Classification of Fusion Experiments:\n")
    
    for category, experiments in classification.items():
        print(f"{category}:")
        for exp in experiments:
            print(f"  - {exp}")
        print() # Add a newline for better spacing

if __name__ == '__main__':
    classify_fusion_experiments()