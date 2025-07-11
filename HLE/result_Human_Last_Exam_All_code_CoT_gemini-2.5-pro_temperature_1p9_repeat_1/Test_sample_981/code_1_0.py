def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    (A) Driving a toroidal current.
    (B) Elongating the flux surfaces and making them rotate poloidally.
    (C) Making the magnetic axis non-planar.

    The interpretation assumes that (A) is a method to achieve the rotational aspect of (B),
    so any device using (A) is also considered to use (B).
    """

    classifications = {
        "Uses only property B": ["Tokamaks", "Reversed Field Pinches"],
        "Uses both B and C": ["LHD", "Wendelstein 7-X"],
        "Uses A, B, and C": ["NCSX"]
    }

    print("Classification of Fusion Experiments by Magnetic Field Twist Mechanism:")
    print("-" * 65)
    print("Property Definitions:")
    print("  (A) Driving a toroidal current")
    print("  (B) Elongating/shaping flux surfaces for poloidal rotation (twist)")
    print("  (C) Having a non-planar magnetic axis")
    print("-" * 65)

    for category, experiments in classifications.items():
        # Format the experiment list for clean printing
        experiment_str = ", ".join(experiments)
        print(f"{category}:")
        print(f"  - {experiment_str}\n")

if __name__ == "__main__":
    classify_fusion_experiments()